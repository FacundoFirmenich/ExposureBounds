"""Exact continuum sign certificates for finite, equally spaced exposure records.

Python standard library only. Coefficients are in ascending power order.
All certificates are conditional mathematical statements, never pilot approvals.
"""
from fractions import Fraction as Q
from math import comb

H0 = (0, 2, 4, 9, 11, 14)
H1 = (0, 2, 5, 7, 12, 14)


def rational(x):
    if isinstance(x, bool) or not isinstance(x, (str, int, Q)):
        raise ValueError('Use integers or exact rational/decimal strings, not floats')
    return Q(x)


def trim(p):
    p = [rational(x) for x in p]
    if not p:
        raise ValueError('Empty polynomial')
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def evaluate(p, x):
    y = Q(0)
    for c in reversed(p):
        y = y*x+c
    return y


def multiply(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def nominal_difference(subdivisions=1, horizon=16):
    if type(subdivisions) is not int or subdivisions < 1:
        raise ValueError('Positive integer subdivisions required')
    if type(horizon) is not int or horizon < 15:
        raise ValueError('Nominal horizon must include all scheduled pulses')
    return [Q(int(i//subdivisions in H1)-int(i//subdivisions in H0))
            for i in range(horizon*subdivisions)]


def filter_difference(difference, kernel, horizon):
    """Exact finite-horizon FIR transport; weights are mathematical model inputs."""
    h = [rational(x) for x in kernel]
    if not h or min(h) < 0 or sum(h) <= 0:
        raise ValueError('A nonzero nonnegative transport kernel is required')
    if type(horizon) is not int or horizon <= 0:
        raise ValueError('Positive integer horizon required')
    return [sum((h[j]*difference[i-j] for j in range(len(h))
                 if 0 <= i-j < len(difference)), Q(0)) for i in range(horizon)]


def coefficients(difference):
    return trim(list(reversed(difference)))


def strip_endpoints(p):
    p = trim(p)
    if p == [0]:
        return p, 0, 0
    zero, one = 0, 0
    while p[0] == 0:
        p = p[1:]
        zero += 1
    while len(p) > 1 and sum(p) == 0:
        # Division by (1-q): b_k = sum_{j<=k} a_j.
        running, quotient = Q(0), []
        for x in p[:-1]:
            running += x
            quotient.append(running)
        assert p[-1] == -quotient[-1]
        p = trim(quotient)
        one += 1
    return p, zero, one


def bernstein(p, a, b):
    """Power expansion under q=a+(b-a)t, then exact power-to-Bernstein."""
    n = len(p)-1
    power = [sum((p[j]*comb(j,k)*a**(j-k)*(b-a)**k
                  for j in range(k,n+1)), Q(0)) for k in range(n+1)]
    return [sum((power[k]*Q(comb(i,k),comb(n,k)) for k in range(i+1)), Q(0))
            for i in range(n+1)]


def certify(p, interval=(Q(0),Q(1)), max_depth=18, max_nodes=8191):
    p = trim(p)
    if len(p) > 129:
        raise ValueError('Degree exceeds the explicit resource limit 128')
    a, b = map(rational, interval)
    if not 0 <= a < b <= 1:
        raise ValueError('Require 0 <= a < b <= 1')
    if type(max_depth) is not int or max_depth < 0 or type(max_nodes) is not int or max_nodes < 1:
        raise ValueError('Invalid certificate resource limits')
    unit_open = a == 0 and b == 1
    if not unit_open and (a == 0 or b == 1):
        raise ValueError('Use the full open unit interval or a closed interior band')
    reduced, z, o = strip_endpoints(p) if unit_open else (p, 0, 0)
    out = {'schema':'living-order-sign-4.9','polynomial':list(map(str,p)),
           'interval':[str(a),str(b)],'domain':'open' if unit_open else 'closed',
           'factors':{'q':z,'one_minus_q':o},'reduced':list(map(str,reduced))}
    if reduced == [0]:
        return dict(out,status='IDENTICALLY_ZERO',leaves=[],lower_bound='0')
    stack, leaves, count = [(a,b,0)], [], 0
    while stack:
        left,right,depth = stack.pop()
        count += 1
        if count > max_nodes:
            return dict(out,status='INCONCLUSIVE',reason='node_budget',nodes=count-1)
        for x in (left,(left+right)/2,right):
            if unit_open and x in (0,1):
                continue
            value = evaluate(p,x)
            if value <= 0:
                return dict(out,status='NEGATIVE_WITNESS' if value < 0 else 'ZERO_WITNESS',
                            witness={'q':str(x),'value':str(value)},nodes=count)
        bs = bernstein(reduced,left,right)
        if min(bs) > 0:
            leaves.append({'interval':[str(left),str(right)],'bernstein':list(map(str,bs))})
        elif depth >= max_depth:
            return dict(out,status='INCONCLUSIVE',reason='depth_budget',nodes=count)
        else:
            mid = (left+right)/2
            stack.extend([(mid,right,depth+1),(left,mid,depth+1)])
    leaves.sort(key=lambda x:Q(x['interval'][0]))
    bound = min(Q(v) for leaf in leaves for v in leaf['bernstein'])
    if z or o:
        bound = Q(0)  # Strict open-interval positivity has no uniform endpoint margin.
    return dict(out,status='POSITIVE',leaves=leaves,lower_bound=str(bound),nodes=count)


def mixture_bound(difference_lower, certificate, step, outside_mass=None):
    """Finite positive spectral measure; mass is normalized only for this bound."""
    step = rational(step)
    if step <= 0:
        raise ValueError('Positive grid step required')
    if certificate['status'] != 'POSITIVE' or certificate['domain'] != 'closed':
        return {'status':'ABSTAIN_NO_POSITIVE_CLOSED_BAND'}
    a = Q(certificate['interval'][0])
    beta = Q(certificate['lower_bound'])
    m = step*a*beta  # (1-q)/(-log q) = integral_0^1 q^u du >= a.
    loss = step*sum((max(Q(0),-rational(x)) for x in difference_lower),Q(0))
    limit = m/(m+loss)
    out = {'band_response_lower':str(m),'global_negative_magnitude_upper':str(loss),
           'outside_mass_strict_threshold':str(limit),
           'normalization':'per unit total nonnegative spectral mass'}
    if outside_mass is None:
        return dict(out,status='ABSTAIN_SPECTRAL_COVERAGE')
    eta = rational(outside_mass)
    if not 0 <= eta <= 1:
        raise ValueError('Outside mass bound must be in [0,1]')
    bound = (1-eta)*m-eta*loss
    return dict(out,status='CONDITIONAL_POSITIVE_MIXTURE' if bound > 0 else 'NOT_CERTIFIED',
                outside_mass_upper=str(eta),mixture_response_lower=str(bound))
