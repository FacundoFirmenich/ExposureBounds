"""Sharp exposure-envelope reduction for positive exponential comparators."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json
from vendor.exact_sign import rational,certify,coefficients,mixture_bound
from vendor.check_certificate import check

V=tuple(map(F,['1','-1/3','-1','-1','0','1','1','1/3','-1']))
def nominal():
    return [F(int(i in {0,2,5,7,12,14})-int(i in {0,2,4,9,11,14})) for i in range(16)]

def worst_profile(epsilon,mean_lower,mean_upper,intra):
    epsilon=rational(epsilon)
    lo,hi,amp=[list(map(rational,a)) for a in (mean_lower,mean_upper,intra)]
    if epsilon<0 or any(len(a)!=16 for a in (lo,hi,amp)):
        raise ValueError('Nonnegative epsilon and sixteen entries required')
    if any(l>u for l,u in zip(lo,hi)) or any(a<0 for a in amp):raise ValueError('Invalid bounds')
    means=nominal()
    for j,v in enumerate(V):means[12-j]-=epsilon*v
    means=[a+b for a,b in zip(means,lo)]
    # Earlier positive, later negative: same extremizer for every positive rate.
    half=[v for m,a in zip(means,amp) for v in (m+a,m-a)]
    return means,half

def analyze(doc):
    if doc.get('schema')!='living-order-continuous-envelope-4.12':raise ValueError('Unknown schema')
    required=['epsilon_coupled_mean','mean_residual_lower','mean_residual_upper','intra_amplitude','Delta_seconds','z_interval']
    missing=[k for k in required if k not in doc]
    if missing:return {'status':'ABSTAIN_MISSING_BOUNDS','missing':missing,'physical_validation':'NOT_PERFORMED'}
    if doc.get('equal_prehistory_declared') is not True or doc.get('envelope_contract_declared') is not True:
        return {'status':'ABSTAIN_INPUT_MODEL','physical_validation':'NOT_PERFORMED'}
    delta=rational(doc['Delta_seconds']);a,b=map(rational,doc['z_interval'])
    if delta<=0 or not 0<a<b<1:raise ValueError('Positive Delta and closed interior z band required')
    if 'outside_spectral_mass_upper' in doc and not 0<=rational(doc['outside_spectral_mass_upper'])<=1:
        raise ValueError('Spectral mass fraction must be in [0,1]')
    means,half=worst_profile(doc['epsilon_coupled_mean'],doc['mean_residual_lower'],doc['mean_residual_upper'],doc['intra_amplitude'])
    p=coefficients(half);global_c=certify(p);band_c=certify(p,(a,b))
    for c in (global_c,band_c):
        if c['status']!='INCONCLUSIVE' and not check(c):raise ValueError('Certificate verification failed')
    if global_c['status']=='POSITIVE':
        spectral={'status':'CONDITIONAL_GLOBAL_POSITIVITY','reason':'All finite positive rates certified; no spectral-band coverage needed.'}
    else:
        spectral=mixture_bound(half,band_c,delta/2,doc.get('outside_spectral_mass_upper'))
    # Bounds on moments of the extra interval-mean error, not on total physical moments.
    lo=list(map(rational,doc['mean_residual_lower']));hi=list(map(rational,doc['mean_residual_upper']))
    moments=[]
    for k in range(3):
        w=[F((i+1)**(k+1)-i**(k+1),k+1)*delta**(k+1) for i in range(16)]
        moments.append([str(sum(x*y for x,y in zip(w,lo))),str(sum(x*y for x,y in zip(w,hi)))])
    return {'schema':'living-order-continuous-result-4.12','input_parameters':doc,
      'worst_interval_means':list(map(str,means)),'worst_half_interval_values':list(map(str,half)),
      'half_step_seconds':str(delta/2),'global_certificate':global_c,'band_certificate':band_c,
      'spectral_conclusion':spectral,'extra_mean_error_moment_bounds':moments,
      'physical_validation':'NOT_PERFORMED',
      'scope':'Exact worst case of the declared measurable exposure family. The half-step witness is not a claim that actual exposure is piecewise constant. Mean residual moments do not include intra-interval moments.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args()
    if args.output.exists():raise ValueError('Refusing to overwrite output')
    raw=args.input.read_bytes();result=analyze(json.loads(raw));result['input_sha256']=hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result.get('spectral_conclusion',result.get('status')))
