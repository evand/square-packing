"""One definition of "same minimum" and of exactsolve verdicts (10-11; REVIEW-2026-10-11.md section 5 item 4).

Verdicts (classify / classify_status), from an exactsolve report or a stored status string:
  certified   exact KKT point with a valid certificate and a second-order proof: 'PD ...', 'rigid ...',
              'local minimum of S, strict modulo N exact flat motions ...'
  weak        valid certificate, second order not proved: 'inconclusive: ...' (negative curvature blocked by weakly
              active contacts) and 'PSD with N zero modes (degenerate ...)'.  Before 10-11 both were counted as
              certified (census_exact.classify); now they are reported separately.
  not-min     'not jammed' (corner-corner MILP descent) or a second-order / first-order 'NOT a local minimum'
  unresolved  anything else (solve failed, timeout, killed, no certificate)

Identity:
  min_key(S)    a minimum = its exact side S (decimal string), quantised to 1e-17.  Sources store S truncated or
                rounded at ~20 decimals; same_min_keys() also returns the two neighbouring keys so a rounding
                difference at the 20th digit cannot split one minimum.  Distinct minima with S closer than 1e-17
                would be an exact tie (same S, different contact graph): not separated here (contact-graph
                tiebreak is a later item).
  side_key(s)   f64 side to 1e-9 (the old basin-table key): a lookup cache for screened sides, NOT a minimum.
                One minimum spans many side keys (unconverged copies sit up to ~1e-3 above it) and two minima can
                share one (10-11 audit: 4 cases at n = 110).
  copy table    NOT identity: 'via descend' rows map a side key to the minimum it descends to (10-20 trap).
  SAME_SIDE     f64 tolerance for "same screened side" (archive dedupe, novelty, chain library): 2e-9.
"""
from __future__ import annotations
from decimal import Decimal, ROUND_HALF_EVEN

CERTIFIED, WEAK, NOT_MIN, UNRESOLVED = 'certified', 'weak', 'not-min', 'unresolved'
RANK = {UNRESOLVED: 0, NOT_MIN: 1, WEAK: 2, CERTIFIED: 3}        # better evidence wins on merge
MIN_OK = ('PD', 'rigid', 'local minimum')
WEAK_PREFIX = ('inconclusive', 'PSD')
SAME_SIDE = 2e-9
_Q = Decimal('1e-17')


def classify_status(status: str | None, has_S: bool = True, cert_valid: bool = True) -> str:
    """Verdict from a (possibly truncated) exactsolve status string."""
    st = (status or '').strip()
    if st.startswith('NOT') or st == 'not jammed':
        return NOT_MIN
    if not (has_S and cert_valid):
        return UNRESOLVED
    if st.startswith(MIN_OK):
        return CERTIFIED
    if st.startswith(WEAK_PREFIX):
        return WEAK
    return UNRESOLVED


def classify(rep: dict) -> str:
    """Verdict from an exactsolve report dict (top-level 'status' set only on failure; 'second_order'.status)."""
    if rep.get('status') == 'not jammed':
        return NOT_MIN
    so = rep.get('second_order') or {}
    st = so.get('status', '') if isinstance(so, dict) else ''
    return classify_status(st, bool(rep.get('S_exact')), bool(rep.get('cert_valid')))


def min_key(S) -> str:
    t = str(S).strip()
    if '/' in t:                                     # rational sides (rational constructions store 'p/q')
        from decimal import localcontext
        num, den = t.split('/')
        with localcontext() as c:
            c.prec = 60
            return str((Decimal(num) / Decimal(den)).quantize(_Q, rounding=ROUND_HALF_EVEN))
    return str(Decimal(t).quantize(_Q, rounding=ROUND_HALF_EVEN))


def same_min_keys(S) -> tuple[str, str, str]:
    k = Decimal(min_key(S))
    return (str(k - _Q), str(k), str(k + _Q))


def side_key(s: float) -> str:
    return f'{float(s):.9f}'


def same_side(a: float, b: float, tol: float = SAME_SIDE) -> bool:
    return abs(a - b) < tol
