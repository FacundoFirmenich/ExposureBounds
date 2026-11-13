"""Reproduce the exact per-arm acquisition example in the supplement."""

from fractions import Fraction
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / 'reproduction'))
from samples_413 import envelope


def main():
    zero, one = Fraction(0), Fraction(1)
    radii = {'A': Fraction(1, 8), 'B': Fraction(1, 16)}
    reports = {}
    for arm, radius in radii.items():
        result = envelope([(zero, one - radius, one + radius)], zero, one)
        assert result['status'] == 'CONSISTENT'
        assert result['pieces'] == [(zero, one, zero, one - radius, zero, one + radius)]
        reports[arm] = {'lower': str(one - radius), 'upper': str(one + radius),
                        'max_input_error': str(radius)}
    # With unit horizon, gain and decay the mean kernel integrates to exp(-1).
    coefficient = sum(radii.values())
    assert coefficient == Fraction(3, 16)
    print(json.dumps({'arms': reports, 'paired_mean_error_bound': '3/(16e)',
                      'coefficient_of_exp_minus_one': str(coefficient)}, indent=2))


if __name__ == '__main__':
    main()
