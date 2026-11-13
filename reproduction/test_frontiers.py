from fractions import Fraction as F
from decimal import Decimal,localcontext
import copy,unittest
from clocks import analyze as clock,verify as vc,shifted
from joint_415 import analyze as previous,tail_certificate
from acquisition_414 import analyze as acquisition
from spectral import analyze as band,verify as vb,exp_minus,response,mixture
from advance import analyze,verify_claims

def m(c,w,l,rho=0):return dict(center=str(c),width=str(w),lower=str(l),upper=str(l),clock_radius=str(rho))
def doc(rows,K='1',shared=None):
    d=dict(schema='living-order-acquisition-4.14',Lipschitz_bound=K,Delta_seconds='1',
        uniform_or_point_sensor_declared=True,simultaneous_bounds_declared=True,equal_prehistory_declared=True,
        clock_model='independent_bounded_offsets',measurements=copy.deepcopy(rows))
    if shared is not None:
        d['clock_model']='shared_offset';d['shared_clock_radius']=str(shared)
        for row in d['measurements']:row.pop('clock_radius')
    return d
def common_case():return doc([m(2,0,0),m(14,0,0),m(8,12,4)],shared=1)
def boundary_case():return doc([m(8,0,0,1),m(8,0,2,1)])
def cover_case():return doc([m(6,0,0),m(10,0,0),m(8,0,3,1)])
def band_case():return doc([m(0,0,3),m(16,0,-1)],K='1/4')

class Tests(unittest.TestCase):
    def test_common_shift_equivalence_closes_previous_abstention(self):
        d=common_case();self.assertEqual(previous(d)['status'],'ABSTAIN_JOINT_UNRESOLVED')
        r=clock(d);self.assertEqual(r['status'],'ORIGINAL_INFEASIBLE');self.assertTrue(vc(d,r))

    def test_nonzero_endpoint_offsets_resolve_boundary_feasibility(self):
        d=boundary_case();self.assertEqual(previous(d,0)['status'],'ABSTAIN_JOINT_UNRESOLVED')
        r=clock(d,max_refinements=0)
        self.assertEqual(r['status'],'ORIGINAL_FEASIBLE');self.assertTrue(vc(d,r))
        self.assertEqual(r['box'],[['-1','-1'],['1','1']])

    def test_cover_excludes_every_clock_value(self):
        d=cover_case();self.assertEqual(previous(d,0)['status'],'ABSTAIN_JOINT_UNRESOLVED')
        r=clock(d,max_refinements=0)
        self.assertEqual(r['status'],'ORIGINAL_INFEASIBLE');self.assertTrue(vc(d,r))
        self.assertEqual(r['boxes_processed'],3)

    def test_missing_cover_branch_rejected(self):
        d=cover_case();r=clock(d,max_refinements=0);r['tree']['children'].pop()
        self.assertFalse(vc(d,r))

    def test_offset_outside_bounds_rejected(self):
        d=boundary_case();r=clock(d,max_refinements=0);r['box'][0]=['-2','-2']
        self.assertFalse(vc(d,r))

    def test_budget_does_not_mean_impossibility(self):
        r=clock(cover_case(),max_boxes=1,max_refinements=0)
        self.assertEqual(r['status'],'ABSTAIN_CLOCK_UNRESOLVED')

    def test_common_equivalence_not_applied_to_independent_offsets(self):
        d=common_case();r=clock(d);other=doc([m(2,0,0,1),m(14,0,0,1),m(8,12,4,1)])
        r['input_parameters']=other
        self.assertFalse(vc(other,r))

    def test_exponential_bounds_against_high_precision_decimal(self):
        with localcontext() as ctx:
            ctx.prec=100
            for x in (F(0),F(1,3),F(1),F(16),F(64)):
                l,u=exp_minus(x);approx=F((-Decimal(x.numerator)/Decimal(x.denominator)).exp())
                self.assertLessEqual(l,approx);self.assertGreaterEqual(u,approx)
                self.assertLessEqual(u-l,F(1,2**80))

    def test_constant_response_analytic_identity(self):
        rows=[['0','16','0','1','0','1']];rate=F(1,8)
        lo,hi=response(rows,rate);e0,e1=exp_minus(2)
        self.assertLessEqual(lo,(1-e0)/rate);self.assertGreaterEqual(hi,(1-e1)/rate)

    def test_band_certifies_despite_negative_tails(self):
        d=band_case();g=acquisition(d)
        self.assertEqual(tail_certificate(g['outer_envelope']['segments'])['status'],'ABSTAIN_TAIL_SIGN')
        r=band(d,['1/100','1/20'])
        self.assertEqual(r['status'],'POSITIVE_ON_BAND_IF_FEASIBLE');self.assertTrue(vb(d,r))

    def test_band_boundary_gap_and_numeric_tampering_rejected(self):
        d=band_case();r=band(d,['1/100','1/20']);r['cells'][0]['interval'][0]='1/99'
        self.assertFalse(vb(d,r))
        r=band(d,['1/100','1/20']);r['cells'][0]['lower']='999'
        self.assertFalse(vb(d,r))

    def test_nonpositive_lower_envelope_does_not_become_original_negative(self):
        r=band(band_case(),['1/2','1'])
        self.assertTrue(r['status'].startswith('ABSTAIN_'))

    def test_resource_domain_abstains(self):
        self.assertEqual(band(band_case(),['5','6'])['status'],'ABSTAIN_EXPONENTIAL_RESOURCE')

    def test_full_claim_verification(self):
        d=band_case();r=analyze(d,['1/100','1/20'],max_refinements=0)
        self.assertEqual(r['clock_analysis']['status'],'ORIGINAL_FEASIBLE');self.assertTrue(verify_claims(d,r))

    def test_band_does_not_imply_unqualified_spectral_mixture(self):
        d=band_case();b=band(d,['1/100','1/20'])
        self.assertEqual(mixture(d,b)['status'],'ABSTAIN_MISSING_SPECTRAL_COVERAGE')
        d['outside_spectral_mass_upper']='1/100'
        self.assertEqual(mixture(d,b)['status'],'ABSTAIN_SPECTRAL_MODEL')

    def test_qualified_spectral_leakage_and_adverse_limit(self):
        d=band_case();d.update(outside_spectral_mass_upper='1/100',nonnegative_normalized_spectrum_declared=True)
        r=analyze(d,['1/100','1/20'],max_refinements=0)
        self.assertEqual(r['spectral_mixture']['status'],'POSITIVE_NORMALIZED_MIXTURE_IF_FEASIBLE');self.assertTrue(verify_claims(d,r))
        d['outside_spectral_mass_upper']='1';b=band(d,['1/100','1/20'])
        self.assertEqual(mixture(d,b)['status'],'ABSTAIN_MIXTURE_SIGN')

    def test_false_mixture_bound_rejected(self):
        d=band_case();d.update(outside_spectral_mass_upper='1/100',nonnegative_normalized_spectrum_declared=True)
        r=analyze(d,['1/100','1/20'],max_refinements=0);r['spectral_mixture']['mixture_lower_bound']='999'
        self.assertFalse(verify_claims(d,r))

if __name__=='__main__':unittest.main(verbosity=2)
