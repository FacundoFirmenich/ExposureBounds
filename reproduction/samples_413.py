"""Exact Lipschitz envelopes from interval-valued samples at exact times."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json
from vendor.exact_sign import rational
from continuous_412 import analyze as certify_exposure,nominal

def cones(samples,K,t):
    return (max(l-K*abs(t-s) for s,l,u in samples),min(u+K*abs(t-s) for s,l,u in samples))

def envelope(samples,K,T):
    samples=sorted(samples)
    if not samples or len(samples)>513 or len({s for s,l,u in samples})!=len(samples):raise ValueError('1..513 distinct sample times required')
    if K<0 or T<=0 or any(not 0<=s<=T or l>u for s,l,u in samples):raise ValueError('Invalid time, interval or regularity bound')
    minimum=max([F(0)]+[max(F(0),l-v)/abs(s-t) for s,l,u in samples for t,a,v in samples if s!=t])
    if K<minimum:return {'status':'INCONSISTENT','minimum_compatible_K':str(minimum)}
    knots=sorted({F(0),T,*[s for s,l,u in samples]});pieces=[]
    for left,right in zip(knots,knots[1:]):
        mid=(left+right)/2
        lower={};upper={}
        for s,l,u in samples:
            sl=-K if mid>=s else K;bl=l-sl*s
            su=K if mid>=s else -K;bu=u-su*s
            lower[sl]=max(lower.get(sl,bl),bl);upper[su]=min(upper.get(su,bu),bu)
        cuts={left,right}
        for lines in (lower,upper):
            if len(lines)==2:
                (a,b),(c,d)=list(lines.items());x=(d-b)/(a-c)
                if left<x<right:cuts.add(x)
        cuts=sorted(cuts)
        for a,b in zip(cuts,cuts[1:]):
            m=(a+b)/2
            sl,bl=max(lower.items(),key=lambda pair:pair[0]*m+pair[1])
            su,bu=min(upper.items(),key=lambda pair:pair[0]*m+pair[1])
            pieces.append((a,b,sl,bl,su,bu))
    sigma=max((u-l)/2 for s,l,u in samples)
    gap_sum=sum((samples[i+1][0]-samples[i][0])**2 for i in range(len(samples)-1))
    width_bound=2*sigma*T+K*(samples[0][0]**2+(T-samples[-1][0])**2+gap_sum/2)
    exact_area=sum((su-sl)*(b*b-a*a)/2+(bu-bl)*(b-a) for a,b,sl,bl,su,bu in pieces)
    assert 0<=exact_area<=width_bound
    return {'status':'CONSISTENT','minimum_compatible_K':str(minimum),'pieces':pieces,
            'integrated_envelope_width':str(exact_area),'resolution_width_bound':str(width_bound)}

def project(pieces,K,delta):
    lo=[];hi=[];amp=[]
    for i in range(16):
        left=i*delta;right=(i+1)*delta;il=F(0);iu=F(0);extlo=[];exthi=[]
        for a,b,sl,bl,su,bu in pieces:
            a=max(a,left);b=min(b,right)
            if a>=b:continue
            il+=sl*(b*b-a*a)/2+bl*(b-a);iu+=su*(b*b-a*a)/2+bu*(b-a)
            extlo.extend((sl*a+bl,sl*b+bl));exthi.extend((su*a+bu,su*b+bu))
        lo.append(il/delta);hi.append(iu/delta)
        amp.append(min(K*delta/2,max(max(exthi)-lo[-1],hi[-1]-min(extlo))))
    return lo,hi,amp

def analyze(doc):
    if doc.get('schema')!='living-order-sampled-exposure-4.13':raise ValueError('Unknown schema')
    if 'Lipschitz_bound' not in doc:return {'status':'ABSTAIN_MISSING_REGULARITY','physical_validation':'NOT_PERFORMED'}
    if doc.get('sample_times_exact_declared') is not True or doc.get('equal_prehistory_declared') is not True:
        return {'status':'ABSTAIN_INPUT_MODEL','physical_validation':'NOT_PERFORMED'}
    K=rational(doc['Lipschitz_bound']);delta=rational(doc['Delta_seconds']);T=16*delta
    points=[tuple(rational(x[k]) for k in ('time_seconds','lower','upper')) for x in doc['samples']]
    env=envelope(points,K,T)
    result={'status':env['status'],'input_parameters':doc,'physical_validation':'NOT_PERFORMED',
            'minimum_compatible_K':env['minimum_compatible_K'],
            'meaning':'K is supplied, not inferred as an upper bound; exact-time deterministic interval contract.'}
    if env['status']!='CONSISTENT':return result
    lo,hi,amp=project(env['pieces'],K,delta);f=nominal()
    bridge={'schema':'living-order-continuous-envelope-4.12','Delta_seconds':str(delta),
      'equal_prehistory_declared':True,'envelope_contract_declared':True,'epsilon_coupled_mean':'0',
      'mean_residual_lower':[str(x-y) for x,y in zip(lo,f)],'mean_residual_upper':[str(x-y) for x,y in zip(hi,f)],
      'intra_amplitude':list(map(str,amp)),'z_interval':doc['z_interval']}
    if 'outside_spectral_mass_upper' in doc:bridge['outside_spectral_mass_upper']=doc['outside_spectral_mass_upper']
    result.update(envelope_segments=[list(map(str,row)) for row in env['pieces']],
      integrated_envelope_width=env['integrated_envelope_width'],resolution_width_bound=env['resolution_width_bound'],
      interval_mean_lower=list(map(str,lo)),interval_mean_upper=list(map(str,hi)),
      intra_amplitude=list(map(str,amp)),bridge_412=bridge,certificate_412=certify_exposure(bridge),
      projection_scope='Pointwise envelopes and mean bounds are sharp; projected product family and amplitude bounds can be conservative.')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args()
    if args.output.exists():raise ValueError('Refusing to overwrite output')
    raw=args.input.read_bytes();r=analyze(json.loads(raw));r['input_sha256']=hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');print(r['status'])
