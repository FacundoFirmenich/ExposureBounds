import json,sys
from advance import verify_claims

if __name__=='__main__':
    with open(sys.argv[1],encoding='utf-8') as f:doc=json.load(f)
    with open(sys.argv[2],encoding='utf-8') as f:r=json.load(f)
    ok=verify_claims(doc,r)
    print('DECISIVE_CLAIMS_VERIFIED_ABSTENTIONS_REMAIN_ABSTENTIONS' if ok else 'NOT_VERIFIED')
    sys.exit(0 if ok else 1)
