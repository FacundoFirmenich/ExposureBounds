"""Clock feasibility plus rate-band certification on the original time domain."""
from pathlib import Path
import argparse,hashlib,json
from clocks import analyze as clocks,verify as verify_clocks
from spectral import analyze as spectral,verify as verify_spectral,mixture
from acquisition_414 import analyze as acquisition
from joint_415 import tail_certificate

def analyze(doc,band=None,max_boxes=31,max_refinements=3,max_cells=64):
    out={'release':'4.16.0','input_parameters':doc,'physical_validation':'NOT_PERFORMED',
         'clock_analysis':clocks(doc,max_boxes,max_refinements),
         'auxiliary_clock_sign_scope':'Nested joint_result tail certificates concern shifted auxiliary contracts only.'}
    g=acquisition(doc)
    if 'outer_envelope' in g:out['original_tail_certificate']=tail_certificate(g['outer_envelope']['segments'])
    if band is not None:
        out['spectral_band']=spectral(doc,band,max_cells)
        out['spectral_mixture']=mixture(doc,out['spectral_band'])
    return out

def verify_claims(doc,r):
    try:
        if r['input_parameters']!=doc:return False
        c=r['clock_analysis']
        if c['status'] in ('ORIGINAL_FEASIBLE','ORIGINAL_INFEASIBLE','INHERITED_GATE'):
            if not verify_clocks(doc,c):return False
        elif c['status']!='ABSTAIN_CLOCK_UNRESOLVED':return False
        g=acquisition(doc)
        if 'outer_envelope' in g and r.get('original_tail_certificate')!=tail_certificate(g['outer_envelope']['segments']):return False
        b=r.get('spectral_band')
        if b and b['status']=='POSITIVE_ON_BAND_IF_FEASIBLE':return verify_spectral(doc,b) and r.get('spectral_mixture')==mixture(doc,b)
        return b is None or b['status'].startswith('ABSTAIN_')
    except (KeyError,ValueError,TypeError):return False

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--band',nargs=2);p.add_argument('--max-boxes',type=int,default=31)
    p.add_argument('--max-refinements',type=int,default=3);p.add_argument('--max-cells',type=int,default=64)
    a=p.parse_args();raw=a.input.read_bytes();d=json.loads(raw)
    r=analyze(d,a.band,a.max_boxes,a.max_refinements,a.max_cells);r['input_sha256']=hashlib.sha256(raw).hexdigest()
    if not verify_claims(d,r):raise ValueError('Internal certificate verification failed')
    with a.output.open('x',encoding='utf-8',newline='\n') as f:json.dump(r,f,indent=2);f.write('\n')
    print(r['clock_analysis']['status'])
