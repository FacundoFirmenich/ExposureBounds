"""Rational enclosures of exponential responses over complete rate intervals."""
from fractions import Fraction as F
from functools import lru_cache
from vendor.exact_sign import rational as Q
from acquisition_414 import analyze as acquisition

@lru_cache(maxsize=4096)
def exp_minus(x):
    x=Q(x)
    if not 0<=x<=64:raise ValueError('Exponential enclosure resource domain: 0..64')
    s=F(1);term=F(1)
    for n in range(513):
        if n+2>x:
            remainder=term*x/(n+1)/(1-x/F(n+2))
            lo,hi=1/(s+remainder),1/s
            if hi-lo<=F(1,2**80):return lo,hi
        term*=x/F(n+1);s+=term
    raise ValueError('Exponential enclosure series budget exhausted')

def multiply(k,interval):
    a,b=interval;return (k*a,k*b) if k>=0 else (k*b,k*a)

def response(segments,rate):
    rate=Q(rate)
    if rate<=0:raise ValueError('Rate must be positive')
    rows=[tuple(map(Q,p)) for p in segments];T=rows[-1][1];lo=hi=F(0)
    for a,b,s,v,_,__ in rows:
        kb=(s*b+v)/rate-s/(rate*rate);ka=(s*a+v)/rate-s/(rate*rate)
        lb,ub=multiply(kb,exp_minus(rate*(T-b)));la,ua=multiply(-ka,exp_minus(rate*(T-a)))
        lo+=lb+la;hi+=ub+ua
    return lo,hi

def derivative_bound(segments):
    rows=[tuple(map(Q,p)) for p in segments];T=rows[-1][1];area=F(0)
    for a,b,s,v,_,__ in rows:
        cuts=[a,b]
        if s and a<-v/s<b:cuts.insert(1,-v/s)
        for l,r in zip(cuts,cuts[1:]):area+=abs(s*(r*r-l*l)/2+v*(r-l))
    return T*area

def record(segments,a,b,B):
    lo,hi=response(segments,(a+b)/2);margin=B*(b-a)/2
    return {'interval':[str(a),str(b)],'lower':str(lo-margin),'upper':str(hi+margin)}

def analyze(doc,band,max_cells=64):
    if type(max_cells)!=int or not 1<=max_cells<=256:raise ValueError('max_cells must be 1..256')
    a,b=map(Q,band)
    if not 0<a<b:raise ValueError('Positive ordered rate band required')
    out={'input_parameters':doc,'band':[str(a),str(b)],'physical_validation':'NOT_PERFORMED','scope':'Conditional on original feasibility; no claim outside this rate band.'}
    g=acquisition(doc)
    if 'outer_envelope' not in g:return dict(out,status='ABSTAIN_ACQUISITION_GATE',gate=g['status'])
    pieces=g['outer_envelope']['segments'];B=derivative_bound(pieces)
    pending=[(a,b)];leaves=[];visited=0
    while pending and visited<max_cells:
        l,r=pending.pop();visited+=1
        try:cell=record(pieces,l,r,B)
        except ValueError as e:return dict(out,status='ABSTAIN_EXPONENTIAL_RESOURCE',reason=str(e))
        if Q(cell['lower'])>0:leaves.append(cell)
        elif Q(cell['upper'])<=0:
            # Only a bound on the lower envelope: NOT a negative original witness.
            return dict(out,status='ABSTAIN_BAND_SIGN',nonpositive_lower_envelope_cell=cell)
        else:
            mid=(l+r)/2;pending.extend(((mid,r),(l,mid)))
    if pending:return dict(out,status='ABSTAIN_BAND_BUDGET',cells_checked=visited)
    return dict(out,status='POSITIVE_ON_BAND_IF_FEASIBLE',derivative_bound=str(B),cells=sorted(leaves,key=lambda p:Q(p['interval'][0])))

def verify(doc,result):
    try:
        if result['input_parameters']!=doc or result['status']!='POSITIVE_ON_BAND_IF_FEASIBLE':return False
        g=acquisition(doc);pieces=g['outer_envelope']['segments'];B=derivative_bound(pieces)
        if result['derivative_bound']!=str(B):return False
        a,b=map(Q,result['band']);cells=result['cells']
        if not 0<a<b or not 1<=len(cells)<=256:return False
        cursor=a
        for cell in cells:
            l,r=map(Q,cell['interval'])
            if l!=cursor or not l<r<=b or cell!=record(pieces,l,r,B) or Q(cell['lower'])<=0:return False
            cursor=r
        return cursor==b
    except (KeyError,ValueError,TypeError,IndexError,ZeroDivisionError):return False

def mixture(doc,band_result):
    """Transfer band positivity only under an explicit normalized spectral contract."""
    if band_result['status']!='POSITIVE_ON_BAND_IF_FEASIBLE':return {'status':'ABSTAIN_NO_POSITIVE_BAND'}
    if 'outside_spectral_mass_upper' not in doc:return {'status':'ABSTAIN_MISSING_SPECTRAL_COVERAGE'}
    if doc.get('nonnegative_normalized_spectrum_declared') is not True:return {'status':'ABSTAIN_SPECTRAL_MODEL'}
    eta=Q(doc['outside_spectral_mass_upper'])
    if not 0<=eta<=1:raise ValueError('Outside spectral mass fraction must be in [0,1]')
    g=acquisition(doc);negative=F(0)
    for row in g['outer_envelope']['segments']:
        a,b,s,v,_,__=map(Q,row);cuts=[a,b]
        if s and a<-v/s<b:cuts.insert(1,-v/s)
        for l,r in zip(cuts,cuts[1:]):negative+=max(F(0),-(s*(r*r-l*l)/2+v*(r-l)))
    margin=min(Q(c['lower']) for c in band_result['cells']);bound=(1-eta)*margin-eta*negative
    return {'status':'POSITIVE_NORMALIZED_MIXTURE_IF_FEASIBLE' if bound>0 else 'ABSTAIN_MIXTURE_SIGN',
            'outside_mass_upper':str(eta),'certified_band_margin':str(margin),'outside_negative_area_bound':str(negative),'mixture_lower_bound':str(bound),
            'scope':'Conditional on declared nonnegative normalized spectrum, coverage bound and original feasibility.'}
