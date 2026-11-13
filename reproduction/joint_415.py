"""Proof-carrying inner/outer Lipschitz feasibility; solver proposals are untrusted."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json
from acquisition_414 import analyze as previous
from vendor.exact_sign import rational as Q

def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))

def model(doc,level,outer):
    if not isinstance(level,int) or isinstance(level,bool) or not 0<=level<=6:
        raise ValueError('Refinement level must be 0..6')
    gate=previous(doc)
    if gate['status'].startswith(('ABSTAIN','INCONSISTENT')):
        raise ValueError('Input must pass the 4.14 necessary gates')
    K=Q(doc['Lipschitz_bound']);T=16*Q(doc['Delta_seconds'])
    rows=[];knots={F(0),T}
    for m in doc['measurements']:
        c,w,l,u=map(Q,(m[k] for k in ('center','width','lower','upper')))
        rho=Q(doc['shared_clock_radius'] if doc['clock_model']=='shared_offset' else m['clock_radius'])
        rows.append((c,w,rho,l,u));knots.update((c-w/2,c+w/2))
    grid=sorted(knots)
    for _ in range(level):grid=sorted(set(grid)|{(a+b)/2 for a,b in zip(grid,grid[1:])})
    if len(grid)>129:raise ValueError('Grid resource limit: 129 nodes')
    n=len(grid);A=[];b=[];errors=[]
    for i,(a,z) in enumerate(zip(grid,grid[1:])):
        v=[F(0)]*n;v[i]=-1;v[i+1]=1
        A.extend((v,[-x for x in v]));b.extend((K*(z-a),K*(z-a)))
    for c,w,rho,l,u in rows:
        v=[F(0)]*n;err=K*rho if outer else F(0)
        if w==0:v[grid.index(c)]=F(1)
        else:
            for i,(a,z) in enumerate(zip(grid,grid[1:])):
                if c-w/2<=a and z<=c+w/2:
                    v[i]+=(z-a)/(2*w);v[i+1]+=(z-a)/(2*w)
                    if outer:err+=K*(z-a)**2/(3*w)
        A.extend((v,[-x for x in v]));b.extend((u+err,-l+err));errors.append(err)
    return grid,A,b,errors

def valid_primal(A,b,x):
    return bool(A) and len(x)==len(A[0]) and all(dot(row,x)<=rhs for row,rhs in zip(A,b))

def valid_dual(A,b,y):
    return len(y)==len(b) and all(v>=0 for v in y) and all(dot(y,col)==0 for col in zip(*A)) and dot(y,b)<0

def solve_equalities(A,b):
    """Exact elimination, setting free variables to zero; None if inconsistent."""
    n=len(A[0]);a=[list(row)+[rhs] for row,rhs in zip(A,b)];piv=[]
    for c in range(n):
        j=next((j for j in range(len(piv),len(a)) if a[j][c]),None)
        if j is None:continue
        r=len(piv);a[r],a[j]=a[j],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for j in range(len(a)):
            if j!=r:
                v=a[j][c];a[j]=[x-v*y for x,y in zip(a[j],a[r])]
        piv.append(c)
    if any(all(x==0 for x in row[:n]) and row[n]!=0 for row in a):return None
    x=[F(0)]*n
    for r,c in enumerate(piv):x[c]=a[r][n]
    return x

def proposal(A,b,dual=False):
    # Floating point is a search device, never a certificate verifier.
    import numpy as np
    from scipy.optimize import linprog
    af=np.array(A,dtype=float);bf=np.array(b,dtype=float)
    options={'time_limit':10.0}
    if not dual:
        sol=linprog(np.zeros(af.shape[1]),A_ub=af,b_ub=bf,bounds=[(None,None)]*af.shape[1],method='highs',options=options)
        if not sol.success:return None,int(sol.status)
        for den in (10**6,10**12):
            x=[F(float(v)).limit_denominator(den) for v in sol.x]
            if valid_primal(A,b,x):return x,int(sol.status)
        # Reconstruct a candidate active face, then recheck EVERY inequality.
        order=np.argsort(np.abs(af@sol.x-bf))
        for count in (af.shape[1],min(len(b),2*af.shape[1])):
            ids=order[:count]
            x=solve_equalities([A[int(i)] for i in ids],[b[int(i)] for i in ids])
            if x is not None and valid_primal(A,b,x):return x,int(sol.status)
        return None,int(sol.status)
    eq=np.vstack((af.T,np.ones(len(b))))
    sol=linprog(bf,A_eq=eq,b_eq=np.r_[np.zeros(af.shape[1]),1.0],bounds=(0,None),method='highs',options=options)
    if not sol.success:return None,int(sol.status)
    for den in (10**6,10**12):
        y=[F(float(v)).limit_denominator(den) for v in sol.x]
        if valid_dual(A,b,y):return y,int(sol.status)
    ids=[i for i,v in enumerate(sol.x) if v>1e-10]
    if ids:
        E=[[A[i][j] for i in ids] for j in range(af.shape[1])]+[[F(1)]*len(ids)]
        v=solve_equalities(E,[F(0)]*af.shape[1]+[F(1)])
        if v is not None:
            y=[F(0)]*len(b)
            for i,a in zip(ids,v):y[i]=a
            if valid_dual(A,b,y):return y,int(sol.status)
    return None,int(sol.status)

def tail_certificate(segments):
    """Exact minima of H(s)=integral_s^T L; segments carry the lower affine lines."""
    right_area=F(0);records=[];minimum=F(0);nonzero=False
    for a,b,s,v,_,__ in reversed([tuple(map(Q,row)) for row in segments]):
        def H(t):return right_area+s*(b*b-t*t)/2+v*(b-t)
        points=[a,b]
        if s and a<-v/s<b:points.append(-v/s)
        values=[H(t) for t in points];minimum=min(minimum,*values)
        nonzero=nonzero or any(x!=0 for x in values)
        records.append({'interval':[str(a),str(b)],'minimum':str(min(values))})
        right_area=H(a)
    return {'status':'ALL_POSITIVE_RATES_IF_ORIGINAL_FEASIBLE' if minimum>=0 and nonzero else 'ABSTAIN_TAIL_SIGN',
            'minimum_tail_area':str(minimum),'total_area':str(right_area),'interval_checks':list(reversed(records))}

def analyze(doc,max_refinements=3):
    if not isinstance(max_refinements,int) or isinstance(max_refinements,bool) or not 0<=max_refinements<=6:
        raise ValueError('max_refinements must be 0..6')
    gate=previous(doc)
    out={'schema':'living-order-joint-result-4.15','physical_validation':'NOT_PERFORMED','input_parameters':doc,'previous_status':gate['status'],'attempts':[]}
    if gate['status'].startswith(('ABSTAIN','INCONSISTENT')):
        return dict(out,status='INHERITED_GATE',gate=gate)
    out['tail_certificate']=tail_certificate(gate['outer_envelope']['segments'])
    for level in range(max_refinements+1):
        try:
            grid,A,b,errors=model(doc,level,False)
        except ValueError as e:
            return dict(out,status='ABSTAIN_RESOURCE_LIMIT',reason=str(e))
        try:x,ps=proposal(A,b)
        except (ImportError,ValueError,OverflowError,FloatingPointError) as e:
            return dict(out,status='ABSTAIN_PROPOSAL_FAILURE',reason=str(e))
        attempt={'level':level,'nodes':len(grid),'primal_solver_status':ps}
        out['attempts'].append(attempt)
        if x is not None:
            cert={'kind':'PRIMAL_PIECEWISE_LINEAR','level':level,'values':list(map(str,x)),'clock_offsets':'ALL_ZERO'}
            out.update(status='ORIGINAL_FEASIBLE',certificate=cert)
            assert verify(doc,out)
            return out
        grid,A,b,errors=model(doc,level,True)
        try:y,ds=proposal(A,b,True)
        except (ImportError,ValueError,OverflowError,FloatingPointError) as e:
            return dict(out,status='ABSTAIN_PROPOSAL_FAILURE',reason=str(e))
        attempt['dual_solver_status']=ds
        attempt['outer_measurement_relaxation']=list(map(str,errors))
        if y is not None:
            out.update(status='ORIGINAL_INFEASIBLE',certificate={'kind':'FARKAS_OUTER','level':level,'multipliers':list(map(str,y)),'contradiction':str(dot(y,b))})
            assert verify(doc,out)
            return out
    return dict(out,status='ABSTAIN_JOINT_UNRESOLVED')

def verify(doc,result):
    """No SciPy dependency: reconstruct constraints from the declared input."""
    try:
        if result['input_parameters']!=doc:return False
        gate=previous(doc)
        if result['status']=='INHERITED_GATE':return gate['status'].startswith(('ABSTAIN','INCONSISTENT')) and result['gate']==gate
        if result.get('tail_certificate')!=tail_certificate(gate['outer_envelope']['segments']):return False
        cert=result['certificate'];outer=cert['kind']=='FARKAS_OUTER'
        grid,A,b,errors=model(doc,cert['level'],outer)
        if outer:
            y=list(map(Q,cert['multipliers']))
            return result['status']=='ORIGINAL_INFEASIBLE' and valid_dual(A,b,y) and cert['contradiction']==str(dot(y,b))
        return cert['kind']=='PRIMAL_PIECEWISE_LINEAR' and result['status']=='ORIGINAL_FEASIBLE' and cert['clock_offsets']=='ALL_ZERO' and valid_primal(A,b,list(map(Q,cert['values'])))
    except (KeyError,ValueError,TypeError,IndexError,ZeroDivisionError):return False

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--max-refinements',type=int,default=3);args=p.parse_args()
    raw=args.input.read_bytes();r=analyze(json.loads(raw),args.max_refinements);r['input_sha256']=hashlib.sha256(raw).hexdigest()
    with args.output.open('x',encoding='utf-8',newline='\n') as f:json.dump(r,f,indent=2);f.write('\n')
    print(r['status'])
