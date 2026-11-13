"""Exact necessary acquisition constraints; no fabricated sensor calibration."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
from vendor.exact_sign import rational as Q
from samples_413 import envelope


def mean_abs(offset, width):
    """E|offset+U| for uniform U on [-width/2,width/2]; Dirac at width=0."""
    d, w = abs(Q(offset)), Q(width)
    if w < 0:
        raise ValueError('Negative width')
    return d if w == 0 or 2*d >= w else w/4+d*d/w


def nullspace(matrix, n):
    """Rational RREF basis of right kernel; no numerical rank tolerance."""
    a = [list(row) for row in matrix]
    pivots = []
    for c in range(n):
        p = next((i for i in range(len(pivots), len(a)) if a[i][c]), None)
        if p is None:
            continue
        r = len(pivots)
        a[r], a[p] = a[p], a[r]
        v = a[r][c]
        a[r] = [x/v for x in a[r]]
        for i in range(len(a)):
            if i != r:
                v = a[i][c]
                a[i] = [x-v*y for x,y in zip(a[i],a[r])]
        pivots.append(c)
    basis = []
    for c in range(n):
        if c in pivots:
            continue
        v = [F(0)]*n
        v[c] = F(1)
        for r,p in enumerate(pivots):
            v[p] = -a[r][c]
        basis.append(v)
    return basis


def functional_relations(rows):
    """Kernel of exact uniform-density and point-mass measurement columns."""
    edges = sorted({x for c,w,r,l,u in rows if w for x in (c-w/2,c+w/2)})
    atoms = sorted({c for c,w,r,l,u in rows if not w})
    matrix = [[F(int(w == 0 and c == t)) for c,w,r,l,u in rows] for t in atoms]
    for a,b in zip(edges,edges[1:]):
        t = (a+b)/2
        matrix.append([1/w if w and c-w/2 < t < c+w/2 else F(0) for c,w,r,l,u in rows])
    return nullspace(matrix,len(rows))


def analyze(doc):
    if doc.get('schema') != 'living-order-acquisition-4.14':
        raise ValueError('Unknown schema')
    out = {'schema':'living-order-acquisition-result-4.14','physical_validation':'NOT_PERFORMED'}
    def finish(status, **fields):
        return dict(out,status=status,**fields)
    if 'Lipschitz_bound' not in doc:
        return finish('ABSTAIN_MISSING_REGULARITY')
    if any(doc.get(k) is not True for k in ('uniform_or_point_sensor_declared','simultaneous_bounds_declared','equal_prehistory_declared')):
        return finish('ABSTAIN_INPUT_MODEL')
    required = ('Delta_seconds','clock_model','measurements')
    if any(k not in doc for k in required):
        return finish('ABSTAIN_MISSING_BOUNDS')
    mode = doc['clock_model']
    if mode not in ('shared_offset','independent_bounded_offsets'):
        return finish('ABSTAIN_CLOCK_MODEL')
    K, delta = Q(doc['Lipschitz_bound']), Q(doc['Delta_seconds'])
    if K < 0 or delta <= 0:
        raise ValueError('K >= 0 and Delta > 0 required')
    if not isinstance(doc['measurements'],list) or not 1 <= len(doc['measurements']) <= 64:
        raise ValueError('1..64 measurements required (implementation resource limit)')
    if mode == 'shared_offset' and 'shared_clock_radius' not in doc:
        return finish('ABSTAIN_MISSING_BOUNDS')
    rows = []
    for m in doc['measurements']:
        keys = ('center','width','lower','upper')
        if any(k not in m for k in keys) or (mode != 'shared_offset' and 'clock_radius' not in m):
            return finish('ABSTAIN_MISSING_BOUNDS')
        if mode == 'shared_offset' and 'clock_radius' in m:
            raise ValueError('Use only shared_clock_radius with shared_offset')
        c,w,l,u = map(Q,(m[k] for k in keys))
        r = Q(doc['shared_clock_radius'] if mode == 'shared_offset' else m['clock_radius'])
        if w < 0 or r < 0 or l > u:
            raise ValueError('Invalid window, clock radius or measurement interval')
        if c-r-w/2 < 0 or c+r+w/2 > 16*delta:
            return finish('ABSTAIN_WINDOW_OUTSIDE_DOMAIN')
        rows.append((c,w,r,l,u))
    out['clock_model'] = mode
    out['input_parameters'] = doc
    out['center_inflation'] = [str(K*mean_abs(r,w)) for c,w,r,l,u in rows]
    pair_checks = []
    for i,(c,w,r,l,u) in enumerate(rows):
        for j,(d,v,s,a,b) in enumerate(rows):
            if i == j:
                continue
            distance = mean_abs(abs(c-d)+(r+s if mode != 'shared_offset' else 0),abs(w-v))
            slack = b+K*distance-l
            if slack < 0:
                return finish('INCONSISTENT_PAIR',witness={'i':i,'j':j,'W1_upper':str(distance),'slack':str(slack)})
            pair_checks.append(slack)
    out['minimum_pair_slack'] = str(min(pair_checks)) if pair_checks else None
    use_relations = mode == 'shared_offset' or all(r == 0 for c,w,r,l,u in rows)
    out['functional_relation_scope'] = 'EXACT_COMMON_TRANSLATION' if use_relations else 'NOT_APPLIED_TO_UNCERTAIN_INDEPENDENT_OFFSETS'
    relations = functional_relations(rows) if use_relations else []
    out['relations_checked'] = []
    for v in relations:
        lo = sum(x*(l if x >= 0 else u) for x,(c,w,r,l,u) in zip(v,rows))
        hi = sum(x*(u if x >= 0 else l) for x,(c,w,r,l,u) in zip(v,rows))
        record = {'coefficients':list(map(str,v)),'interval':[str(lo),str(hi)]}
        out['relations_checked'].append(record)
        if not lo <= 0 <= hi:
            return finish('INCONSISTENT_LINEAR_FUNCTIONALS',witness=record)
    merged = {}
    for c,w,r,l,u in rows:
        radius = K*mean_abs(r,w)
        a,b = l-radius,u+radius
        if c in merged:
            a,b = max(a,merged[c][0]),min(b,merged[c][1])
        if a > b:
            return finish('INCONSISTENT_CENTER_BOUNDS',witness={'center':str(c)})
        merged[c] = (a,b)
    env = envelope([(c,a,b) for c,(a,b) in merged.items()],K,16*delta)
    if env['status'] != 'CONSISTENT':
        return finish('INCONSISTENT_CENTER_BOUNDS',witness=env)
    pieces = env.pop('pieces')
    out['outer_envelope'] = dict(env,segments=[list(map(str,p)) for p in pieces],scope='OUTER_RELAXATION_NOT_SHARP_FOR_ORIGINAL_AVERAGES')
    lower_min = min(s*t+v for a,b,s,v,_,__ in pieces for t in (a,b))
    lower_area = sum(s*(b*b-a*a)/2+v*(b-a) for a,b,s,v,_,__ in pieces)
    out['direct_sign'] = {'lower_min':str(lower_min),'lower_area':str(lower_area),'status':'CONDITIONAL_ALL_POSITIVE_RATES_IF_FEASIBLE' if lower_min >= 0 and lower_area > 0 else 'ABSTAIN_SIGN','scope':'Sufficient positive-kernel bound; failure is not a negative witness.'}
    low,high = max(l for c,w,r,l,u in rows),min(u for c,w,r,l,u in rows)
    if low <= high:
        return finish('FEASIBLE_CONSTANT_WITNESS',constant_witness=str(low),clock_offset_witness='0')
    return finish('OUTER_FEASIBLE_ORIGINAL_UNRESOLVED')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args()
    raw=args.input.read_bytes();result=analyze(json.loads(raw))
    result['input_sha256']=hashlib.sha256(raw).hexdigest()
    with args.output.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(result['status'])
