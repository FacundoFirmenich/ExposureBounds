"""Compare verified subdivisions for the release's formal spectral example."""
import json
from fractions import Fraction as F

from acquisition_414 import analyze as acquisition
from spectral import derivative_bound, exp_minus, record
from test_frontiers import band_case


def local_bound(segments, rate):
    """Rational upper enclosure using exponential and lag maxima per segment."""
    horizon = F(segments[-1][1])
    bound = F(0)
    for raw in segments:
        left, right, slope, intercept, _, _ = map(F, raw)
        cuts = [left, right]
        if slope and left < -intercept / slope < right:
            cuts.insert(1, -intercept / slope)
        for a, b in zip(cuts, cuts[1:]):
            area = abs(slope * (b*b - a*a) / 2 + intercept * (b - a))
            _, upper_weight = exp_minus(rate * (horizon - b))
            bound += (horizon - a) * upper_weight * area
    return min(bound, derivative_bound(segments))


def compare(segments, left, right, bound, max_cells=256):
    pending = [(left, right)]
    leaves = []
    visited = 0
    while pending and visited < max_cells:
        a, b = pending.pop()
        visited += 1
        cell = record(segments, a, b, bound(segments, a))
        if F(cell['lower']) > 0:
            leaves.append(cell)
        elif F(cell['upper']) <= 0:
            return {'status': 'ABSTAIN_ENVELOPE', 'checked': visited}
        else:
            mid = (a + b) / 2
            pending.extend(((mid, b), (a, mid)))
    if pending:
        return {'status': 'ABSTAIN_BUDGET', 'checked': visited}
    leaves.sort(key=lambda cell: F(cell['interval'][0]))
    cursor = left
    for leaf in leaves:
        a, b = map(F, leaf['interval'])
        assert a == cursor and a < b <= right
        assert leaf == record(segments, a, b, bound(segments, a))
        cursor = b
    assert cursor == right
    return {'status': 'CERTIFIED', 'checked': visited, 'leaves': len(leaves),
            'margin': float(min(F(cell['lower']) for cell in leaves))}


def main():
    segments = acquisition(band_case())['outer_envelope']['segments']
    left, right = F(1, 100), F(1, 20)
    baseline = lambda pieces, _rate: derivative_bound(pieces)
    result = {'adaptive': {
        'global': compare(segments, left, right, baseline),
        'local': compare(segments, left, right, local_bound),
    }}
    for cells in (1, 4, 16):
        grid = [left + (right-left) * F(k, cells) for k in range(cells + 1)]
        result[f'fixed_{cells}'] = {
            label: float(min(F(record(segments, a, b, bound(segments, a))['lower'])
                             for a, b in zip(grid, grid[1:])))
            for label, bound in [('global', baseline), ('local', local_bound)]
        }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
