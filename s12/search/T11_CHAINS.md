# T11_CHAINS: every tight wall-to-wall chain of the verified `T = 11` packing, and which hypothesis of the chain lemma fails there  (2026-09-21)

*Diagnostic, one hour.  Input: `runs/bandcut_scan_exact110.json` (the exact `n = 110`, `T = 11` certificate of
`search/BANDCUT_SCAN.md` §2, `delta = +2.748946598e-04`).  Script `search/t11_chains.py`; output `runs/t11_chains_110.json`.
Everything here is **[measured]** on one configuration except the inequality in §1, which is **[proved]** (three lines).*

## 0. Verdict

1. **Every tight wall-to-wall chain runs through the band.**  There is no wall-to-wall path of axis-parallel squares at all
   (the band cuts every row and column, as `ARCH_TLEDGER.md` said of Cleemann).  At slack `< 1e-6` there are exactly 2 tight
   x-chains and 3 tight y-chains (up to which axis-parallel row they use), all of length exactly `T = 11`, each with 2–4 tilted
   squares at one end; the rest of the chain is an axis-parallel row.
2. **The hypothesis that fails is rise control, not the chain's existence and not the mixed tilts.**  On each tight chain the
   rise-free budget `L` is strongly negative (`-0.62` to `-0.96`; down to `-3.2` on near-tight chains) and the rise credit `C`
   is positive and slightly larger.  The chain descends by `2.4–4.7` units over the container — the `R* ~ 1.1` of
   `BANDCUT_K.md` §1.2 is exceeded four-fold — and nothing transverse pins it, because the transverse neighbours of a band
   square are other band squares whose links are themselves tilted staircases.
3. **The margin is a difference of two `O(1)` numbers per chain.**  Per tilted link the credit is `tan(tau) * |Dperp|` with
   `|Dperp| ~ 0.5` (the band steps half a unit sideways per link) and `tau ~ 25–29°`, i.e. `+0.24–0.28`, against a length
   cost `m/cos(tau) - 1 ~ 0.13–0.16`: net **`+0.10–0.13` per band link**.  The two interfaces cost `-0.16` (the tilt-change
   link, `(W(D)-1)/2` of `bandcut-cost.md` §1) and `-0.18` (a tilted end square against the wall, `(u-1)/2`).  Three band
   links at `+0.116` pay `+0.35` against `-0.34`: the whole `+2.7e-4` is that residue, divided by `D ~ 12.5`.
4. **So the `T`-dependent quantity, in chain language, is: how many full-credit band links a wall-to-wall chain can contain
   between its two interfaces.**  At `T = 11` a chain of eleven has room for 3–4 band links *and* seven axis-parallel squares
   and two interfaces.  At `T = 4` a chain of four that contains three band links is the band itself, with no axis-parallel
   part and its interfaces on the walls — the `(4,1)` stack of `BANDCUT_SCAN.md` §5 — which is a margin-`0` configuration
   still cut by an axis-parallel chain.  **[heuristic]**  This is the same conclusion as `bandcut-cost.md` §2.2 ("the band
   must be thick") with the number attached: a band link is worth `~+0.12`, an interface `~-0.17`, and a chain has to carry
   `>= 3` band links to break even.

## 1. The inequality used (exact, per-link normals)

For a path `Q_1..Q_k` whose link `s` is separated along an x-type unit normal `n_s = (cos tau_s, sin tau_s)` (an edge normal
of either square of the link, `cos tau_s > 0`), with `m_s = 1/2 + W(theta_{s+1} - theta_s)/2` and `p = W(theta)/2`:

    delta * D  <=  L + C,     D = 2 + sum_s 1/cos tau_s,
    L = T - p_1 - p_k - sum_s m_s / cos tau_s          (rise-free budget; `0` for `T` axis-parallel squares in a row),
    C = sum_s tan(tau_s) * Dy_s                        (rise credit; the only term that can be positive).

Divide each link row by `cos tau_s`, add the two wall rows, sum the x-components.  **[proved]**  At the configuration itself
`(L + C)/D - delta = [wall slacks + sum_s (gap_s - delta)/cos tau_s]/D >= 0`; that quantity is the **slack** of the path, and
a slack-`0` path is a tight chain.  MT (`bandcut-cost.md` §1) is the special case of a common normal; the per-link form
telescopes for any normals because it projects on the *global* axis, at the price of the `tan(tau_s) Dy_s` terms — which is
where the rise enters.  Per link the net contribution is `1 - m_s/cos tau_s + tan(tau_s) Dy_s`, `0` for an axis-parallel link.

## 2. The tight chains

`python3 search/t11_chains.py runs/bandcut_scan_exact110.json --top 400 --cap 1e-6 --detail 30`.  65 near-axis (`< 1°`) and
45 tilted squares.  Distinct tight chains (the axis-parallel row they continue into is free; `Dperp` is the sideways step of
each link, `net` the per-link contribution of §1):

| axis | tilts along the chain | link normals | `Dperp` on tilted links | `net` per link | `L` | `C` | rise |
|---|---|---|---|---|---|---|---|
| x | `-29.6 -29.5 -29.5 -25.4` then 7 × `0` | `-29.5 -29.5 -29.5 -25.4`, then `0` | `-.48 -.49 -.50 -.27` | `+.122 +.130 +.096 -.162` | `-0.960` | `+0.963` | `-4.75` |
| x | `-24.9 -24.9 -24.9` then 8 × `0` | `-24.9 × 3`, then `0` | `-.50 -.43 -.48` | `+.131 +.097 -.061` | | | |
| y | 8 × `0` then `+23.9 +23.9 +23.9` | `0`, then `-23.9 × 3` | `-.47 -.45 -.47` | `-.060 +.107 +.116` | `-0.616` | `+0.620` | `-2.40` |
| y | 7 × `0` then `+27.3 +28.1 +28.1 +28.1` | `0`, then `-27.3 -28.1 × 3` | `-.30 -.48 -.47 -.47` | `-.169 +.112 +.118 +.118` | | | |
| y | 9 × `0` then `+28.3 +28.3` | `0`, then `-28.3 × 2` | `-.76 -.45` | `+.073 +.108` | | | |

Every tight chain has `(L+C)/D = delta` to `1e-8` (the rational rounding of the certificate).  Near-tight census (`--cap`):

| slack `<=` | x-chains | y-chains | lengths seen | tilted squares per chain |
|---|---|---|---|---|
| `1e-6` | 2 | 3 | 11 | 2–4 |
| `1e-3` | 10 | 38 | 10–12 | 2–10 |
| `3e-2` | 48 | 102 | 9–12 | 1–10 |

No wall-to-wall path of axis-parallel squares exists at any slack.  Over all chains of slack `< 1e-3`: `L in [-3.2, -0.4]`,
`C in [0.4, 3.2]`, `|Dperp|` on tilted links mean `0.50` (min `0.001`, max `1.12`).

## 3. Reading against the candidate lemma

The candidate of `notes/proof-architecture.md` §0a item 14 — "a tight wall-to-wall chain of `T` on which MT gives `<= 0`" —
is satisfied in its first half (tight chains of `T` exist, five of them) and fails in the second: MT needs a common normal
and these chains change normal by `24–29°` at the interface, and the per-link form of §1, which has no such restriction,
gives exactly `delta > 0` on them because it is tight.  The content of a chain certificate is therefore entirely in
**bounding `C`**, i.e. the `T`-shaped H-lemma of `tasks/t3-chain`: at `T = 11` the rise of the band part is `~0.5` per link
and unpinned.  What a `T = 4` proof has to use is that a chain of four cannot afford three band links plus two interfaces
(item 4 above) — a counting statement in the chain length, which is where `4` enters as a number.  **[heuristic]**

Not done: the same listing on the `T = 12` certificate (`runs/bandcut_scan_exact132.json`; same script, one command);
transverse-chain search through the band squares to quantify "unpinned".
