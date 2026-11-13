"""Verified coverage of clock-offset boxes and common-translation reduction."""
from copy import deepcopy
from itertools import product
from collections import deque
from vendor.exact_sign import rational as Q
from acquisition_414 import analyze as acquisition
from joint_415 import analyze as joint,verify as verify_joint

def reject(r):
    return r['status']=='ORIGINAL_INFEASIBLE' or (r['status']=='INHERITED_GATE' and r['gate']['status'].startswith('INCONSISTENT'))

def root_box(doc):
    radii=[Q(doc['shared_clock_radius'])] if doc['clock_model']=='shared_offset' else [Q(m['clock_radius']) for m in doc['measurements']]
    return [(-r,r) for r in radii]

def shifted(doc,box):
    d=deepcopy(doc)
    for i,m in enumerate(d['measurements']):
        a,b=box[0 if d['clock_model']=='shared_offset' else i]
        m['center']=str(Q(m['center'])+(a+b)/2)
        if d['clock_model']!='shared_offset':m['clock_radius']=str((b-a)/2)
    if d['clock_model']=='shared_offset':d['shared_clock_radius']=str((box[0][1]-box[0][0])/2)
    return d

def encoded(box):return [[str(a),str(b)] for a,b in box]

def analyze(doc,max_boxes=31,max_refinements=3,probe_corners=True):
    if type(max_boxes)!=int or not 1<=max_boxes<=127:raise ValueError('max_boxes must be 1..127')
    if type(max_refinements)!=int or not 0<=max_refinements<=6:raise ValueError('max_refinements must be 0..6')
    gate=acquisition(doc)
    out={'input_parameters':doc,'physical_validation':'NOT_PERFORMED','schema':'living-order-clock-result-4.16'}
    if gate['status'].startswith(('ABSTAIN','INCONSISTENT')):
        return dict(out,status='INHERITED_GATE',gate=gate)
    bounds=root_box(doc)
    if doc['clock_model']=='shared_offset':
        # Existence is translation invariant for this exact acquisition contract.
        zero=[(Q(0),Q(0))];r=joint(shifted(doc,zero),max_refinements)
        status='ORIGINAL_FEASIBLE' if r['status']=='ORIGINAL_FEASIBLE' else 'ORIGINAL_INFEASIBLE' if reject(r) else 'ABSTAIN_CLOCK_UNRESOLVED'
        return dict(out,status=status,method='COMMON_TRANSLATION_EQUIVALENCE',joint_result=r)
    probes=0
    if probe_corners:
        active=[i for i,(a,b) in enumerate(bounds) if a!=b]
        if len(active)<=4:
            for bits in product((0,1),repeat=len(active)):
                point=[(Q(0),Q(0)) for _ in bounds]
                for i,bit in zip(active,bits):point[i]=(bounds[i][bit],bounds[i][bit])
                r=joint(shifted(doc,point),max_refinements);probes+=1
                if r['status']=='ORIGINAL_FEASIBLE':
                    return dict(out,status='ORIGINAL_FEASIBLE',method='EXPLICIT_OFFSET_WITNESS',box=encoded(point),joint_result=r,corner_probes=probes)
    tree={};queue=deque([(tree,bounds)]);processed=0
    while queue and processed<max_boxes:
        node,box=queue.popleft();r=joint(shifted(doc,box),max_refinements);processed+=1
        if r['status']=='ORIGINAL_FEASIBLE':
            return dict(out,status='ORIGINAL_FEASIBLE',method='EXPLICIT_OFFSET_WITNESS',box=encoded(box),joint_result=r,boxes_processed=processed,corner_probes=probes)
        if reject(r):node.update(kind='REJECTED',joint_result=r);continue
        widths=[b-a for a,b in box];i=max(range(len(box)),key=lambda i:widths[i])
        if widths[i]==0:
            node.update(kind='UNRESOLVED',joint_result=r);continue
        a,b=box[i];mid=(a+b)/2;left=list(box);right=list(box)
        left[i]=(a,mid);right[i]=(mid,b)
        node.update(kind='SPLIT',axis=i,children=[{},{}])
        queue.extend(((node['children'][0],left),(node['children'][1],right)))
    for node,box in queue:node.update(kind='UNVISITED')
    def closed(node):return node.get('kind')=='REJECTED' or (node.get('kind')=='SPLIT' and all(closed(c) for c in node['children']))
    return dict(out,status='ORIGINAL_INFEASIBLE' if closed(tree) else 'ABSTAIN_CLOCK_UNRESOLVED',method='EXHAUSTIVE_BOX_COVER',tree=tree,boxes_processed=processed,corner_probes=probes)

def verify(doc,result):
    try:
        if result['input_parameters']!=doc:return False
        gate=acquisition(doc)
        if result['status']=='INHERITED_GATE':return gate['status'].startswith(('ABSTAIN','INCONSISTENT')) and result['gate']==gate
        if gate['status'].startswith(('ABSTAIN','INCONSISTENT')):return False
        bounds=root_box(doc);method=result['method']
        if method=='COMMON_TRANSLATION_EQUIVALENCE':
            if doc['clock_model']!='shared_offset':return False
            d=shifted(doc,[(Q(0),Q(0))]);r=result['joint_result']
            return verify_joint(d,r) and ((result['status']=='ORIGINAL_FEASIBLE' and r['status']=='ORIGINAL_FEASIBLE') or (result['status']=='ORIGINAL_INFEASIBLE' and reject(r)))
        if method=='EXPLICIT_OFFSET_WITNESS':
            box=[tuple(map(Q,p)) for p in result['box']]
            if len(box)!=len(bounds) or any(not a<=c<=d<=b for (a,b),(c,d) in zip(bounds,box)):return False
            r=result['joint_result']
            return result['status']=='ORIGINAL_FEASIBLE' and r['status']=='ORIGINAL_FEASIBLE' and verify_joint(shifted(doc,box),r)
        if method!='EXHAUSTIVE_BOX_COVER' or result['status']!='ORIGINAL_INFEASIBLE':return False
        count=[0]
        def check(node,box):
            count[0]+=1
            if count[0]>255:return False
            if node['kind']=='REJECTED':
                r=node['joint_result'];return reject(r) and verify_joint(shifted(doc,box),r)
            if node['kind']!='SPLIT' or len(node['children'])!=2:return False
            i=node['axis']
            if type(i)!=int or not 0<=i<len(box):return False
            a,b=box[i]
            if a>=b:return False
            mid=(a+b)/2;left=list(box);right=list(box);left[i]=(a,mid);right[i]=(mid,b)
            return check(node['children'][0],left) and check(node['children'][1],right)
        return check(result['tree'],bounds)
    except (KeyError,ValueError,TypeError,IndexError,RecursionError):return False
