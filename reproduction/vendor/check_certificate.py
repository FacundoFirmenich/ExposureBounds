"""Independent certificate checker: expands Bernstein bases into global powers.

Does not import the prover, its interval transformation, or its subdivision code.
"""
from fractions import Fraction as F
from math import comb
import argparse,json
from pathlib import Path

def mul(a,b):
    r=[F(0)]*(len(a)+len(b)-1)
    for i in range(len(a)):
        for j in range(len(b)): r[i+j]+=a[i]*b[j]
    return r

def power(p,n):
    r=[F(1)]
    for _ in range(n): r=mul(r,p)
    return r

def clean(p):
    while len(p)>1 and p[-1]==0: p.pop()
    return p

def check(c):
    if c.get('schema')!='living-order-sign-4.9': raise ValueError('Unknown schema')
    original=list(map(F,c['polynomial']));reduced=list(map(F,c['reduced']))
    z,o=c['factors']['q'],c['factors']['one_minus_q']
    if type(z) is not int or type(o) is not int or min(z,o)<0: raise ValueError('Bad factors')
    reconstructed=mul(mul(power([0,1],z),power([1,-1],o)),reduced)
    if clean(reconstructed)!=clean(original[:]): raise ValueError('Factorization mismatch')
    a,b=map(F,c['interval'])
    if not 0<=a<b<=1: raise ValueError('Domain mismatch')
    if c['domain']=='open':
        if (a,b)!=(0,1): raise ValueError('Unsupported open domain')
    elif c['domain']=='closed':
        if not 0<a<b<1 or z or o: raise ValueError('Invalid closed band')
    else: raise ValueError('Invalid domain kind')
    status=c['status']
    if status in ('NEGATIVE_WITNESS','ZERO_WITNESS'):
        x=F(c['witness']['q']);v=sum((v*x**i for i,v in enumerate(original)),F(0))
        valid=(a<x<b) if c['domain']=='open' else (a<=x<=b)
        if not valid or str(v)!=c['witness']['value']: raise ValueError('Invalid witness')
        if (status=='NEGATIVE_WITNESS' and v>=0) or (status=='ZERO_WITNESS' and v!=0):
            raise ValueError('Invalid witness sign')
        return True
    if status=='IDENTICALLY_ZERO':
        if any(original): raise ValueError('Not identically zero')
        return True
    if status=='INCONCLUSIVE':
        return False  # No sign theorem is being certified.
    if status!='POSITIVE': raise ValueError('Unknown status')
    cursor=a;minimum=None
    if not c['leaves']: raise ValueError('No coverage')
    for leaf in c['leaves']:
        l,r=map(F,leaf['interval']);bs=list(map(F,leaf['bernstein']));n=len(bs)-1
        if l!=cursor or not l<r<=b or n!=len(reduced)-1 or min(bs)<=0:
            raise ValueError('Invalid coverage or positivity')
        polynomial=[F(0)]*(n+1)
        for k,v in enumerate(bs):
            term=mul(power([-l,1],k),power([r,-1],n-k))
            scale=v*comb(n,k)/(r-l)**n
            for j,x in enumerate(term): polynomial[j]+=scale*x
        if clean(polynomial)!=clean(reduced[:]): raise ValueError('Bernstein identity mismatch')
        minimum=min(bs) if minimum is None else min(minimum,min(bs))
        cursor=r
    if cursor!=b: raise ValueError('Incomplete interval coverage')
    expected=F(0) if z or o else minimum
    if F(c['lower_bound'])!=expected: raise ValueError('Invalid uniform lower bound')
    return True

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('certificate',type=Path)
    args=p.parse_args();ok=check(json.loads(args.certificate.read_text(encoding='utf-8')))
    print('VERIFIED' if ok else 'NO_SIGN_CERTIFICATE');raise SystemExit(0 if ok else 2)
