#!/usr/bin/env python3
"""Convert a plain-text certificate to a self-describing ``points.json`` and back.

The plain-text format (``certificates/FORMAT.md``) is deliberately bare: whitespace-
separated integers, nothing else.  This script publishes the same data in the JSON
schema used by Mira's ``17squares`` repo (``points.json``), extended with ``weights`` /
``weight_denominator`` and with provenance, so that the certificate is readable by
either family of checkers.  The JSON is a *companion*, never the source of truth: the
``.txt`` is what ``verify/`` reads and what ``certificates/SHA256SUMS`` pins.

Both directions are exact (integers only; the ``decimal`` fields are informational and
are ignored on import) and deterministic, so::

    txt -> json -> txt   reproduces the shipped .txt byte for byte, and
    txt -> json          reproduces the shipped .json byte for byte.

Usage
-----
    python3 search/export_points.py export   CERT.txt  [-o CERT.json] [--n 12]
    python3 search/export_points.py import   CERT.json [-o CERT.txt]
    python3 search/export_points.py --roundtrip CERT.txt CERT.json

``--roundtrip`` is the check ``verify.sh`` runs: it asserts that the JSON regenerates
the shipped text file byte for byte, that the text file regenerates the shipped JSON
byte for byte, that the SHA-256 recorded in the JSON is the text file's, and that the
recorded total weight is the sum of the weights and is below ``n``.  Only the standard
library is used.
"""

import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

SCHEMA = "s12-points-json/1"
DEFAULT_N = 12

CONVENTION = {
    "squares": "closed",
    "boundary_points_count": True,
    "claim": (
        "every closed unit square contained in the closed container [0, s]^2, "
        "at every centre and every angle, contains points of total weight >= 1"
    ),
    "argument": (
        "concentric shrink: a packing of n unit squares in a container of side s' < s, "
        "rescaled by s/s' > 1, gives n squares each strictly containing a concentric "
        "closed unit square; those have pairwise disjoint interiors, so no point is "
        "counted twice and n <= total_weight < n is a contradiction"
    ),
    "note": (
        "this differs from the open-square / pigeonhole convention of the unweighted "
        "family (Fort, Mira): a certificate for one convention does not verify under "
        "the other's checker"
    ),
    "reference": "certificates/FORMAT.md",
}


# --------------------------------------------------------------------------- text side

def read_txt(text):
    """Parse certificate text. Returns (s: Fraction, D, W, rows: list[(x, y, w)])."""
    toks = text.split()
    it = iter(toks)
    s_num, s_den = int(next(it)), int(next(it))
    D, W, m = int(next(it)), int(next(it)), int(next(it))
    rows = [(int(next(it)), int(next(it)), int(next(it))) for _ in range(m)]
    rest = list(it)
    if rest:
        raise ValueError(f"{len(rest)} trailing tokens after {m} points")
    return Fraction(s_num, s_den), D, W, rows


def write_txt(s, D, W, rows):
    """Canonical text form, identical to search/scale_to_critical.py:write_cert."""
    out = [f"{s.numerator} {s.denominator}\n{D}\n{W}\n{len(rows)}\n"]
    out += [f"{x} {y} {w}\n" for x, y, w in rows]
    return "".join(out)


# --------------------------------------------------------------------------- json side

def rational(fr):
    """Exact rational as {numerator, denominator, decimal}; decimal is informational."""
    fr = Fraction(fr)
    return {
        "numerator": fr.numerator,
        "denominator": fr.denominator,
        "decimal": float(fr),
    }


def txt_to_dict(txt_bytes, txt_name, n):
    s, D, W, rows = read_txt(txt_bytes.decode("ascii"))
    total = Fraction(sum(w for _, _, w in rows), W)
    if not total < n:
        raise ValueError(f"total weight {total} is not below n = {n}; not a certificate")
    bound = f"{s.numerator}/{s.denominator}"
    return {
        "schema": SCHEMA,
        "problem": f"packing {n} unit squares in a square",
        "theorem": f"s({n}) >= {bound}",
        "bound": bound,
        "n": n,
        "container_side": rational(s),
        "point_denominator": D,
        "weight_denominator": W,
        "points": [[x, y] for x, y, _ in rows],
        "weights": [w for _, _, w in rows],
        "total_weight": rational(total),
        "convention": CONVENTION,
        "source": {
            "file": txt_name,
            "sha256": hashlib.sha256(txt_bytes).hexdigest(),
            "format": "certificates/FORMAT.md",
        },
    }


def dict_to_txt(d):
    """Rebuild the canonical text certificate from the exact integer fields only."""
    s = Fraction(d["container_side"]["numerator"], d["container_side"]["denominator"])
    D, W = int(d["point_denominator"]), int(d["weight_denominator"])
    pts, wts = d["points"], d["weights"]
    if len(pts) != len(wts):
        raise ValueError(f"{len(pts)} points but {len(wts)} weights")
    rows = [(int(x), int(y), int(w)) for (x, y), w in zip(pts, wts)]
    return write_txt(s, D, W, rows)


def dumps(obj, indent=0):
    """Deterministic pretty-printer: 2-space indent, one point / one weight per line."""
    pad = " " * indent
    if isinstance(obj, dict):
        items = [f'{pad}  {json.dumps(k)}: {dumps(v, indent + 2)}' for k, v in obj.items()]
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    if isinstance(obj, list):
        if len(obj) <= 4 and all(not isinstance(x, (list, dict)) for x in obj):
            return json.dumps(obj)
        items = [f"{pad}  {dumps(v, indent + 2)}" for v in obj]
        return "[\n" + ",\n".join(items) + "\n" + pad + "]"
    if isinstance(obj, bool) or obj is None or isinstance(obj, (int, str)):
        return json.dumps(obj)
    if isinstance(obj, float):
        return repr(obj)
    raise TypeError(type(obj))


def to_json_bytes(d):
    return (dumps(d) + "\n").encode("ascii")


# --------------------------------------------------------------------------- commands

def cmd_export(args):
    txt = Path(args.txt)
    d = txt_to_dict(txt.read_bytes(), txt.name, args.n)
    out = Path(args.output) if args.output else txt.with_suffix(".json")
    out.write_bytes(to_json_bytes(d))
    print(f"wrote {out}  ({len(d['points'])} points, total weight "
          f"{d['total_weight']['numerator']}/{d['total_weight']['denominator']})")


def cmd_import(args):
    src = Path(args.json)
    d = json.loads(src.read_bytes())
    out = Path(args.output) if args.output else src.with_suffix(".txt")
    out.write_bytes(dict_to_txt(d).encode("ascii"))
    print(f"wrote {out}")


def roundtrip(txt_path, json_path):
    txt_path, json_path = Path(txt_path), Path(json_path)
    T, J = txt_path.read_bytes(), json_path.read_bytes()
    d = json.loads(J)
    fails = []

    # 1. json -> txt reproduces the shipped text file byte for byte.
    if dict_to_txt(d).encode("ascii") != T:
        fails.append("json -> txt does not reproduce the shipped .txt byte for byte")

    # 2. txt -> json reproduces the shipped json byte for byte (n taken from the json).
    if to_json_bytes(txt_to_dict(T, txt_path.name, int(d["n"]))) != J:
        fails.append("txt -> json does not reproduce the shipped .json byte for byte")

    # 3. provenance: the recorded hash and file name are the text file's.
    h = hashlib.sha256(T).hexdigest()
    if d["source"]["sha256"] != h:
        fails.append(f"source.sha256 {d['source']['sha256']} != sha256(.txt) {h}")
    if d["source"]["file"] != txt_path.name:
        fails.append(f"source.file {d['source']['file']!r} != {txt_path.name!r}")

    # 4. arithmetic: total_weight is the sum of the weights, and is below n.
    total = Fraction(sum(int(w) for w in d["weights"]), int(d["weight_denominator"]))
    rec = Fraction(d["total_weight"]["numerator"], d["total_weight"]["denominator"])
    if total != rec:
        fails.append(f"total_weight {rec} != sum of weights {total}")
    if not total < int(d["n"]):
        fails.append(f"total weight {total} is not below n = {d['n']}")
    if d["bound"] != f"{d['container_side']['numerator']}/{d['container_side']['denominator']}":
        fails.append("bound does not match container_side")

    tag = f"{txt_path.name} <-> {json_path.name}"
    if fails:
        print(f"ROUNDTRIP FAILED  {tag}")
        for f in fails:
            print(f"  - {f}")
        return False
    print(f"roundtrip OK  {tag}  ({len(d['points'])} points, total weight {rec}, "
          f"sha256 {h[:16]}...)")
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--roundtrip", nargs=2, metavar=("CERT.txt", "CERT.json"),
                    help="assert the two files are exact images of each other")
    sub = ap.add_subparsers(dest="cmd")
    ex = sub.add_parser("export", help="txt -> json")
    ex.add_argument("txt")
    ex.add_argument("-o", "--output")
    ex.add_argument("--n", type=int, default=DEFAULT_N,
                    help=f"number of squares the certificate rules out (default {DEFAULT_N})")
    ex.set_defaults(func=cmd_export)
    im = sub.add_parser("import", help="json -> txt")
    im.add_argument("json")
    im.add_argument("-o", "--output")
    im.set_defaults(func=cmd_import)
    args = ap.parse_args(argv)
    if args.roundtrip:
        return 0 if roundtrip(*args.roundtrip) else 1
    if not args.cmd:
        ap.print_help()
        return 2
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
