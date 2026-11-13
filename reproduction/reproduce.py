"""Formal clock and spectral proof objects; not experimental observations."""
import json
from pathlib import Path
from test_frontiers import common_case,boundary_case,cover_case,band_case
from advance import analyze,verify_claims

def run():
    covered=band_case();covered.update(outside_spectral_mass_upper='1/100',nonnegative_normalized_spectrum_declared=True)
    cases={
       'common_offset_infeasible':analyze(common_case()),
       'nonzero_boundary_offsets_feasible':analyze(boundary_case(),max_refinements=0),
       'all_offsets_excluded':analyze(cover_case(),max_refinements=0),
       'band_positive_with_negative_tails':analyze(band_case(),['1/100','1/20'],max_refinements=0),
       'unresolved_budget_preserved':analyze(cover_case(),max_boxes=1,max_refinements=0),
       'qualified_spectral_mixture':analyze(covered,['1/100','1/20'],max_refinements=0),
    }
    for r in cases.values():assert verify_claims(r['input_parameters'],r)
    assert [r['clock_analysis']['status'] for r in cases.values()]==['ORIGINAL_INFEASIBLE','ORIGINAL_FEASIBLE','ORIGINAL_INFEASIBLE','ORIGINAL_FEASIBLE','ABSTAIN_CLOCK_UNRESOLVED','ORIGINAL_FEASIBLE']
    assert cases['band_positive_with_negative_tails']['spectral_band']['status']=='POSITIVE_ON_BAND_IF_FEASIBLE'
    return {'release':'4.16.0','evidence_kind':'FORMAL_PROOF_OBJECTS_NOT_MEASUREMENTS','cases':cases}

if __name__=='__main__':
    r=run();p=Path('results');p.mkdir(exist_ok=True)
    (p/'FORMAL_RESULTS.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v['clock_analysis']['status'] for k,v in r['cases'].items()},indent=2))
