"""jlevy register witnesses with algebraic-number-field coordinates (format since 2026-11-05, register #913)."""
import math
from pk import packing as P

W = """softschema:
  contract: packing.squares:Witness/v2
witness:
  coordinates: {angle_unit: not-applicable, axes: x-right-y-up, origin: lower-left}
  n: 2
  representation: center-basis
  scalar:
    kind: algebraic-number-field
    minimal_polynomial: ['1', '0', '-2']
    isolating_interval: ['1', '2']
  side: ['1/2', '1']
  squares:
  - basis: [['1', '0'], ['0', '0']]
    center: [['1/2', '0'], ['1/2', '0']]
  - basis: [['1', '0'], ['0', '1']]
    center: [['0', '1/2'], ['1', '1/4']]
"""


def test_field_center_basis():
    p = P.parse(W)
    assert abs(p.s - (0.5 + math.sqrt(2))) < 1e-15
    assert abs(p.sq[1][0] - math.sqrt(2) / 2) < 1e-15 and abs(p.sq[1][1] - (1 + math.sqrt(2) / 4)) < 1e-15
    assert abs(p.sq[0][2]) < 1e-12 and abs(p.sq[1][2] - math.degrees(math.atan2(math.sqrt(2), 1))) < 1e-9


def test_field_val_higher_degree():
    v = P.field_val({'minimal_polynomial': ['1', '0', '0', '-2'], 'isolating_interval': ['1', '2']})
    a = 2 ** (1 / 3)
    assert abs(v(['1', '1', '1']) - (1 + a + a * a)) < 1e-14 and v('3/4') == 0.75
