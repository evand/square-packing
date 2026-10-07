# friedman-lit — literature check (2026-09-30; subagent report, saved by coordinator, condensed)

Labels: [P] read in primary source, [S] secondhand, [I] inference.  Roth–Vaughan (JCTA 24 (1978) 170–186) itself was NOT
reachable (403s); all statements of it are secondhand.

## 1. Roth–Vaughan
* Abstract [P-ish, PSU record]: "in packing a square of side n + 1/2 with unit squares, the wasted space always has area ≫ n^{1/2}."
* Friedman DS7 (1998 pdf and 2009 html agree) [S]: "if s(s − ⌊s⌋) > 1/6, then W(s) ≥ 10^{−100} √(s · |s − ⌊s + .5⌋|)".
  So **W(α) ≥ 10⁻¹⁰⁰ √(α·‖α‖)**, ‖α‖ = distance to nearest integer.  Constant 10⁻¹⁰⁰ (10⁻¹⁰ is RV's tilt threshold).
  Wikipedia: Ω((a·|a − round a|)^{1/2}); Bui 2025: W(x) ∉ o(x^{1/2}).  McClenagan 2026 has a garbled restatement.
* **Consequence [I]: RV does NOT give s(k²−c) = k for k ≥ k0(c).**  Side α = k − δ holding k² − c squares forces kδ ≲ c/2,
  and RV only gives W ≥ 10⁻¹⁰⁰ √(kδ) ≲ 10⁻¹⁰⁰ √c ≪ c.  The bound vanishes as α → k⁻.  (QUADRANT.md §7(iii) is wrong.)
  Only if the paper's bound were √(α(α−⌊α⌋)) would it work; no source states that.  Check Theorem 1 of the paper if ever in hand.
* Method [S, via Bui arXiv 2504.09489 which restates RV's lemmas]: "good" squares (tilt ≤ 10⁻¹⁰) vs bad; fundamental lemma
  "difference in angle leads to wasted space" (a disk of radius 2 around a square near a differently-tilted neighbour/wall
  contains waste ≥ c·θ); counting along vertical lines/strips.  Geometric rigidity, not an LP/weighting argument; no obvious
  fractional analogue [I].  Bui Thm 1: W(x) = Θ(W*(x)), W* = packings with good squares only.

## 2. Friedman DS7 conjectures [P]
* 1998 v1: "**Conjecture 1.** If s(n²−k)=n, then s((n+1)²−k)=n+1. That is, if omitting k squares from an n × n square does
  not admit a smaller packing, then the same will be true for omitting k squares from any larger perfect square packing.
  This is true of all the best known packings."  "**Conjecture 2.** W(s) = O(s^{1/2})."
* 2000 and 2009 versions: Conjecture 1 only (Conj. 2 dropped).  Evidence given: only "true of all the best known packings".
* Also: "It was conjectured that s(n²−n)=n whenever n is small. The smallest known counterexample ... Cleemann, s(17²−17)<17."

## 3. Work bearing on Conj. 1 / c*(k) = max{c : s(k²−c)=k}
* Nagamochi EJC 12 (2005) R37: s(n²−2)=s(n²−1)=n ∀n ⇒ c*(k) ≥ 2.
* Bentz EJC 17 (2010) R126 (s(13)=4, s(46)=7); arXiv 1606.03746 (s(22)=5, s(33)=6), "strongly suggest s(m²−3)=m".
* Arslanov, Mustafin, Shangitbayev, EJC 28(4) (2021) P4.22, doi 10.37236/8586 [P pp.1–2]: s(n²−n) < n for all n ≥ 12
  (⇒ c*(k) ≤ k−1 for k ≥ 12); also s(18²−17)<18, s(17²−16)<17, s(16²−15)<16.  Strip lemma for rectangles:
  δ((Rx,Ry),m) ≤ δ((Rx+1,Ry), m+Ry−1) ("adding a unit strip costs one square") — the ascent direction.
  They note Erdős–Graham ⇒ s(n² − O(n^{7/11})) < n.
* **No published result that c*(k) → ∞.**  No paper/MO/blog attacking Conjecture 1 found (Ellsworth, MathWorld,
  Wikipedia, Friedman mirror checked).

## 4. Waste W(s) upper bounds
1975 Erdős–Graham O(s^{7/11}); 1978 Montgomery (unpubl.) O(s^{(3−√3)/2+ε}); 2009 Chung–Graham O(s^{(3+√2)/7} log s);
2016 Wang–Dong–Li arXiv 1603.02368 O(s^{5/8}); 2020 Chung–Graham DCG 64 claimed O(s^{3/5}) — error (Arslanov–Bui note, DCG
2025, doi 10.1007/s00454-025-00767-w); 2025 Bui arXiv 2508.04603 O(s^{0.6}); 2026 McClenagan arXiv 2602.01484 O(s^{3/5}).
Lower: only Roth–Vaughan.
URLs: DS7 1998 https://www.combinatorics.org/files/Surveys/ds7/ds7v1-1998.pdf ; 2009
https://www.combinatorics.org/files/Surveys/ds7/ds7v5-2009/ds7-2009.html ; Ellsworth https://kingbird.myphotos.cc/packing/squares.html
