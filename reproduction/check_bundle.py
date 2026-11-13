import hashlib,json,sys
from pathlib import Path
from advance import verify_claims

def check(actual_path,reference_path):
    a=Path(actual_path).read_bytes();b=Path(reference_path).read_bytes();actual=json.loads(a);reference=json.loads(b)
    assert actual['release']==reference['release']=='4.16.0'
    assert actual['evidence_kind']==reference['evidence_kind']=='FORMAL_PROOF_OBJECTS_NOT_MEASUREMENTS'
    assert actual['cases'].keys()==reference['cases'].keys()
    for k,r in actual['cases'].items():
        old=reference['cases'][k];d=r['input_parameters']
        assert d==old['input_parameters'] and verify_claims(d,r) and verify_claims(d,old)
        assert r['clock_analysis']['status']==old['clock_analysis']['status']
        assert r.get('original_tail_certificate')==old.get('original_tail_certificate')
        assert r.get('spectral_band')==old.get('spectral_band')
        assert r.get('spectral_mixture')==old.get('spectral_mixture')
    return {'all_decisive_claims_verified':True,'same_inputs_and_conclusions':True,'full_bytes_equal':a==b,
      'actual_sha256':hashlib.sha256(a).hexdigest(),'reference_sha256':hashlib.sha256(b).hexdigest(),'physical_validation':'NOT_PERFORMED'}

if __name__=='__main__':print(json.dumps(check(sys.argv[1],sys.argv[2]),indent=2))
