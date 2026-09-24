# Deformations of Global Isolation in the Lonely Runner Conjecture

**T. Palelei.** v9.0.

---

## Abstract

Runners with distinct constant speeds move round a circular track of length 1, and a runner is
**lonely** at a moment when every other runner is at least `1/n` away. The Lonely Runner Conjecture
says every runner is lonely at some moment. When the speeds are equally spaced, all `n` runners are
lonely at once, evenly spread round the track, and this happens exactly `φ(n)` times before the
pattern repeats (Palelei 2026a). We call this **global isolation** and take it as the ideal. Any
other speeds are the ideal plus deviations: we give the exact displacement each deviation causes,
and the exact condition under which a moment of global isolation survives it. When it does not, each
runner's isolation is shifted, and we show where it goes. For any runner of any arrangement, with no
hypothesis beyond distinct speeds, the runner is lonely at some moment exactly when it is lonely at a
moment of global isolation of the equally spaced arrangement whose step is one of its own speed
differences `e` — a moment `(j + 1/n)/e`, called an **opening**. A deterministic construction,
moved only by the runners in the way, reaches such an opening exactly when the runner is lonely. For
whole-number speeds, whether an opening is lonely is a check of `n − 1` remainders. In these terms
the conjecture says the isolation is only ever shifted, never destroyed, and a proof is a finite list
of conditions on the speed differences, each with a proof that a lonely opening exists under it. We
give seven such conditions at general order, and show how they compose: differences not divisible
by `n` of small enough total cost can be removed, passing the question to a smaller set at the same
order. The results are formalised in Lean 4.

---

## 1. Introduction

Let `v₁, …, vₙ` be distinct real speeds, each runner moving on a circle of unit circumference at its
own speed, and let `‖x‖` denote the distance from `x` to the nearest integer. Runner `i` is **lonely**
at time `t` if `‖(vⱼ − vᵢ)t‖ ≥ 1/n` for every `j ≠ i`. The Lonely Runner Conjecture, due to Wills and
independently to Cusick, asserts that every runner is lonely at some moment. It is proved for
`n ≤ 7` by arguments particular to each order, the last two being Bohman, Holzman and Kleitman (2001)
at six runners and Barajas and Serra (2008) at seven. It has been verified by computer at eight
runners (Rosenfeld 2025) and, by extensions of that method, up to thirteen (Sungkawichai and
Trakulthongchai 2026). Nothing in this paper depends on those results, or on any result outside the
framework developed here.

**Global isolation is the ideal.** If the runners' speeds are equally spaced, say `1, 2, …, n`, there
are moments when all `n` runners are lonely at once: they stand evenly spread around the track, `1/n`
apart. We call this **global isolation**. It happens exactly `φ(n)` times before the whole pattern
repeats (Theorem 4); with eight runners, four times, at `t = 1/8, 3/8, 5/8` and `7/8`. Any other set
of speeds is the ideal plus deviations, and the deviations displace each runner from its evenly
spread position by an exact amount (Theorem 5). Theorem 6 says exactly which deviations keep a moment
of global isolation: every runner must be displaced by a whole number of steps of `1/n`, and no two
runners may land on the same point. When that fails, the runners no longer line up together, and each
runner's isolation is shifted to a moment of its own. We show where it goes: a runner is lonely
exactly when it is lonely at a moment of global isolation of an equally spaced arrangement whose step
is one of its own speed differences (§§5–6). The conjecture, in these terms, is that the isolation is
only ever shifted: no arrangement of speeds can shift it out of existence.

We make no claim about *which* lonely moment is found — not the first, not the loneliest. The claim
is only that the isolation of the ideal is shifted by a change of speeds, deterministically, to a
moment of a known form, and is not destroyed.

**How the shifted isolation is found.** Fix one runner and call it the subject. Another runner is *in
the way* while it is closer to the subject than `1/n`, and stops being in the way the moment it is
exactly `1/n` away. So each other runner is in the way for a short stretch every time it passes the
subject, and each stretch ends with that runner exactly `1/n` away. We prove that if the subject is
ever lonely, it is lonely at one of those end moments (Theorem 9). Those end moments are the openings,
and each is a moment of global isolation of the equally spaced arrangement whose step is the
difference between the subject's speed and that runner's (§6). For whole-number speeds, checking an
opening is a short calculation with remainders (Theorem 10). The result is a chain of equivalences,
each proved in both directions for every arrangement:

```
    loneliness   ⟺   lonely opening   ⟺   certificate   ⟺   cascade halts
```

Here a **certificate** is a pair `(e, j)` naming a lonely opening, and the **cascade** is a
construction that starts when the nearest runner first leaves the subject and moves only when some
runner is in the way (§7).

**What is proved.** Because the chain is an equivalence, a condition on the speed differences under
which a lonely opening exists is a proof of the conjecture for every arrangement meeting it — the
isolation is shown to be shifted and not destroyed there, with nothing further owed. We prove seven
such conditions at every order (§10), and show how they compose (§11): when the differences not
divisible by `n` have small enough total cost, removing them passes the question to a smaller set at
the same order, and any condition that settles the smaller set settles the original.

**What remains** is to close the list: conditions of this kind that jointly cover every set of
differences (§14). When the list closes, the isolation is shown never to be destroyed.

**Organisation.** Layer 1, §§2–8, establishes the characterisation. Layer 2, §§9–14, gives the
conditions under which a lonely opening exists and describes what they leave. §15 records what holds
at general `n`, and §16 the formal verification.

---

# LAYER 1 — THE FRAMEWORK

---

## 2. Loneliness is a closed condition, and its boundary

Let `v₁, …, vₙ` be distinct real speeds and fix a **subject** `k`. Runner `i` sits at circle distance
`‖(vᵢ − v_k)t‖` from the subject at time `t`. Write the subject's **disparities** for its positive
speed differences `|vᵢ − v_k|`; there are at most `n − 1` of them.

> **Definition.** The subject is **lonely at `t`** when `‖d·t‖ ≥ 1/n` for every disparity `d`.

The condition is **closed**: a runner at distance exactly `1/n` leaves the subject alone. Everything
below runs on this.

At any instant the `n` runners cut the circle into `n` cyclic gaps summing to 1.

> **Theorem 1.** Runner `i` is lonely if and only if both gaps flanking it are at least `1/n`.

*Proof.* Order the runners around the circle, and write `g_j` for the gap between consecutive
runners `j` and `j+1`, indices mod `n`.

Suppose both gaps flanking `i` are at least `1/n`, and let `m ≠ i`. The circle is cut by `i` and `m`
into two arcs whose lengths sum to 1, and the circle distance is the shorter of them. One arc runs
from `i` in the direction of increasing index and contains the whole of `g_i`; the other runs in the
direction of decreasing index and contains the whole of `g_{i−1}`. Both arcs are therefore at least
`1/n`, hence so is the shorter.

Conversely, suppose runner `i` is lonely, and take `m = i+1`. One of the two arcs between them is
exactly `g_i`, and the circle distance is the smaller of the two arcs, so it is at most `g_i`.
Loneliness gives `1/n ≤ g_i`. Taking `m = i−1` gives `1/n ≤ g_{i−1}`. ∎

> **Theorem 2.** All `n` runners are simultaneously lonely if and only if every gap is exactly
> `1/n` — the regular `n`-gon.

*Proof.* If every gap is exactly `1/n` then every runner's two flanking gaps are `1/n`, so every
runner is lonely by Theorem 1. Conversely, suppose all `n` are lonely. Each gap flanks the two
runners at its ends, so Theorem 1 gives `g_j ≥ 1/n` for every `j`. The gaps sum to 1, so the slacks
`g_j − 1/n` are nonnegative and sum to 0, and each is zero. ∎

> **Corollary 3.** The number of simultaneously lonely runners lies in `{0, …, n−2} ∪ {n}`. It is
> never exactly `n − 1`.

*Proof.* Suppose exactly `n − 1` runners are lonely, all but runner `j`. Every gap has two endpoints,
at most one of them `j`, so every gap flanks a lonely runner and is at least `1/n` by Theorem 1. The
gaps sum to 1, so every gap is exactly `1/n`, both gaps flanking `j` are `1/n`, and `j` is lonely too.
∎

Theorem 2 says simultaneous isolation is rigid. Corollary 3 is a conservation statement: the gaps sum
to 1, so isolation cannot be removed from one runner while the others keep theirs simultaneously.

### 2.1 The boundary discipline

Because the condition is closed, the set of times at which a given runner obstructs is a union of
**open** intervals,

```
    ( (N − 1/n)/d , (N + 1/n)/d ),     N ∈ ℤ,
```

one about each of that runner's relative laps. We call these its **in-the-way zones**; the endpoints
belong to the complement.

**Endpoints are clear.** At `t = (N + 1/n)/d` the runner `d` stands at distance exactly `1/n`. Every
opening, and every instant the cascade of §7 visits, is such an endpoint for some runner and some
lap — which is what makes the certificate of §6 a pair of integers rather than a real number.

**A cell is decided by an interior point.** Partition the positive axis by all the zone endpoints of
all the disparities. Within one cell every runner's status is constant, so any interior point decides
it. The endpoints bounding a cell may both be clear while the cell's interior is not: the two bounding
zones can belong to different runners and overlap across the cell. Deciding a cell therefore requires
an interior point.

This is the **boundary set method** of Palelei (2026a, §2.2), whose Lemmas 2.1–2.2 establish that a
runner's status can change only at these endpoints: evaluate open cells and boundary points
separately, never merge touching zones, and decide a cell by an interior point. The example of §7.1
shows it operating.

### 2.2 What the method uses, and what a proof is

The arguments below use exact arithmetic and the geometry of lap zones. They use no measure, no
density, no equidistribution and no covering, and each exclusion has a reason in the object.

**No measure or density argument.** At the arithmetic progression — which Theorem 2 identifies as
the rigid configuration any general argument must cover — the instants at which all `n` runners are
lonely are the `φ(n)` points per period of Theorem 4, a set of measure zero. There are no
neighbourhoods there for an averaging argument to act on.

**No covering as the form of a conclusion.** The construction forms no complement and exhibits no
set to be shown empty; it runs to a named instant and produces the pair `(e, j)`. Counting arguments
— pigeonhole, residue counts — are admissible as tools, and several proofs below use them.

**No bound on when, and no count of steps.** The conjecture asks for loneliness at *some* positive
time. A bound on when the instant falls, or on how many steps of the cascade reach it, is a stronger
obligation than the conjecture carries.

**No search.** The cascade is fixed by the speeds: the runners in the way at `t` determine the move,
and nothing is ranked or guessed (§7).

**No global clock, and no pinned runner.** Each runner receives its own instant, read in its own
differences `v_j − v_i`; no runner is placed at rest.

These choices fix the form of a proof. A proof in this framework is a finite list of cases, each a
condition on the disparity set together with a proof that a certificate exists under it, the list
closing when the conditions are jointly exhaustive. Each case is discharged by establishing
existence — by pigeonhole, induction, contradiction or construction.

---

## 3. Global isolation, and how it deforms

Theorem 2 says that simultaneous isolation is rigid. There is one family of arrangements at which
the configuration is reached, its instants can be written down, and there are `φ(n)` of them in each
period; every other arrangement is a deformation of it. The `φ(n)` law is prior work (Palelei 2026a,
Theorem 3.1), restated with its proof so the paper is self-contained, and no part of it is claimed as
new. What is developed here is the deformation: an exact criterion for when a moment of global isolation survives
(§3.1), and, where none does, the location of the isolation that remains (§§4–8).

> **Theorem 4** (the `φ(n)` law; Palelei 2026a, Theorem 3.1)**.** Let the speeds be
> `v₀, v₀ + d, …, v₀ + (n−1)d`. All `n` runners are simultaneously lonely precisely at the instants
>
> ```
>     t = a/(n d),      a ≥ 1,      gcd(a, n) = 1,
> ```
>
> and at no others. There are `φ(n)` such instants in each period `1/d`.

*Proof* (Palelei 2026a; restated). Work relative to runner `0`. Writing `u = d t`, the `n` relative
positions are `{ i u mod 1 : 0 ≤ i < n }`. By Theorem 2 all `n` runners are lonely exactly when these
are the vertices of a regular `n`-gon containing `0`, that is

```
    { i u mod 1 : 0 ≤ i < n }  =  { 0, 1/n, 2/n, …, (n−1)/n }.
```

If they are, runner `1` stands at a vertex, so `u = a/n` with `a ≥ 1`. As `i` runs over `ℤ/n` the
image of `i ↦ i a` is the subgroup generated by `gcd(a, n)`, of size `n / gcd(a, n)`, which is all of
`ℤ/n` exactly when `gcd(a, n) = 1`. Conversely, if `gcd(a, n) = 1` then `i ↦ i a` is a bijection of
`ℤ/n` and the positions are precisely the `n` vertices. ∎

### 3.1 Deformation, and why the problem goes runner by runner

Palelei (2026a, Proposition 4.1) shows that the configuration of Theorem 4 has zero stability
radius. What follows sharpens that to an exact criterion.

Let the speeds be arbitrary. Write them against a progression, `v = A + e`, with `A_k = v₀ + k d` and
deviations `e_k`, normalising `e₀ = 0`.

> **Theorem 5.** At `t = a/(n d)` with `gcd(a, n) = 1`, the progression part places the `n` runners
> one to each `n`-th; the deviation displaces runner `k` from its vertex by exactly `e_k · a / d` of an
> `n`-th.

*Proof.* At `t = a/(n d)` runner `k` sits at

```
    (v₀ + k d + e_k) · a/(n d)  =  v₀ a/(n d)  +  (k a)/n  +  (e_k a / d) · (1/n).
```

The first term is a common rotation. The second places runner `k` at the vertex `k a / n`, and these
are all `n` vertices since `gcd(a,n) = 1`. The third is the displacement, of size `e_k a / d` in units
of `1/n`. ∎

> **Theorem 6.** The moment of global isolation survives the deviation if and only if every displacement is a
> whole number of `n`-ths and `k ↦ k a + e_k a / d (mod n)` is a bijection of `ℤ/n`.

*Proof.* Write `w_k = e_k a / d`, so that runner `k` sits at `(k a + w_k)/n` up to the common
rotation.

Suppose each `w_k` is an integer and `k ↦ k a + w_k` is a bijection of `ℤ/n`. For `k ≠ k'` the two
runners are separated by `((k a + w_k) − (k' a + w_{k'}))/n`, whose numerator is an integer not
divisible by `n`, so their circle distance is at least `1/n`. Every runner is lonely.

Conversely, suppose all `n` runners are lonely at `t`. By Theorem 2 the positions are a translate of
the grid `{0, 1/n, …, (n−1)/n}`, so every `k a + w_k` differs from an integer by a quantity
independent of `k`; since `w₀ = 0` that quantity is `0`, and every `w_k` is an integer. The `n`
positions are distinct grid points, so `k ↦ k a + w_k (mod n)` is injective on `ℤ/n`, hence a
bijection. ∎

So simultaneity requires the displacements to be integral, which forces commensurable speeds, and
requires the resulting map of the vertices to remain a permutation. An irrational deviation breaks the
global moment outright; a rational one breaks it as soon as two runners are sent to the same vertex.

**Where the isolation goes.** When the moment of global isolation fails, only the simultaneity is
lost: the runners no longer line up together, and each runner's isolation is shifted to a moment of
its own. By Corollary 3 this cannot happen to one runner alone — if one runner is not lonely at a
moment, at least one other is not lonely there either. Layer 1 shows where each runner's isolation is
shifted to. Layer 2 proves, for every condition on its list, that it is not destroyed.

---

## 4. One runner, one scale

Normalise by the nearest runner. Write `d_min` for the least disparity, `s = d / d_min ≥ 1` for the
**ratio** of a disparity `d`, and `u = n·d_min·t` for the local coordinate. The **reading** of `d` at
time `t` is `s·u = n·d·t`, and `d` is clear exactly when its reading lies at distance at least `1`
from `n·ℤ`.

> **Theorem 7 (the earliest possible lonely instant).** If the subject is lonely at `t > 0` then
> `t ≥ 1/(n·d_min)`; equivalently `u ≥ 1`.

*Proof.* If `0 < t < 1/(n·d_min)` then `0 < d_min·t < 1/n`, so `‖d_min·t‖ < 1/n`. ∎

We write `t₀ = 1/(n·d_min)` and call it **the bound**: below it the nearest runner has not yet left
the subject. It is where the construction of §7 begins. The whole problem is now visible in one line:
**find `u ≥ 1` at which every reading lies at distance at least 1 from `n·ℤ`.**

> **Theorem 8 (the ratio condition).** If `d_max ≤ (n−1)·d_min` then the subject is lonely at the
> bound `t₀ = 1/(n·d_min)`.

*Proof.* At `u = 1` every reading is its ratio `s ∈ [1, n−1]`, at distance at least 1 from `n·ℤ`. ∎

The constant `n − 1` is one below the first multiple of `n`. Theorem 8 returns in Layer 2 as rule R1.

---

## 5. Localisation: every lonely instant is represented by an opening

> **Definition.** For a disparity `e` and an integer `j ≥ 0`, the `j`-th **opening** of the flank
> `e` is the instant `(j + 1/n)/e`.

At an opening the runner `e` sits at distance exactly `1/n` — on the boundary, hence clear.

> **Theorem 9 (Localisation).** The subject is lonely at some positive instant if and only if some
> flank has a lonely opening.

*Proof.* An opening is a positive instant, so one direction is immediate.

For the other let `t > 0` be lonely. For a disparity `e` put `T_e = (⌊e·t⌋ + 1/n)/e`, the opening of
`e` at its last completed lap before `t`. Since `‖e·t‖ ≥ 1/n`, the fractional part of `e·t` is at
least `1/n`, which says exactly `T_e ≤ t`. Let `T = max_e T_e`, attained at some `e*`, so that `T` is
an opening of `e*` and `T ≤ t`.

`T` is lonely. Let `d` be any disparity. From `T ≥ T_d` we get `d·T ≥ ⌊d·t⌋ + 1/n`, and from `T ≤ t`
we get `d·T ≤ d·t < ⌊d·t⌋ + 1`. So the fractional part of `d·T` is `d·T − ⌊d·t⌋ ∈ [1/n, {d·t}]`, and
`{d·t} ≤ 1 − 1/n` because `t` is lonely. So `‖d·T‖ ≥ 1/n`. ∎

The proof produces, from any lonely `t`, a lonely opening no later than `t`. So every lonely instant
is represented by a lonely opening, and everything after this works with the countable family of
openings and loses nothing by doing so.

---

## 6. The certificate

> **Definition.** A **certificate** for a subject is a pair `(e, j)` — a flank and an opening index —
> such that the opening `(j + 1/n)/e` is lonely.

By Theorem 9, the subject is lonely at some instant exactly when it has a certificate. For
whole-number disparities the condition is arithmetic.

> **Theorem 10.** Let the disparities be positive integers and `e` one of them. The opening `j` of
> the flank `e` is lonely if and only if
>
> ```
>     d · (n j + 1)  mod  n e   ∈   [ e , (n−1) e ]        for every disparity d.
> ```
>
> The condition is periodic in `j` with period `e`, so the `e` indices `0 ≤ j < e` decide the
> question outright.

*Proof.* At the opening `(j + 1/n)/e` the disparity `d` stands at `d(n j + 1)/(n e)`. Distance at
least `1/n` from the integers says, after multiplying by `n e`, that `x = d(n j + 1) mod n e`
satisfies `x ≥ e` and `n e − x ≥ e`. For periodicity,
`d(n(j + e) + 1) = d(n j + 1) + d·n e`. ∎

Written on the floor rather than the residue, with `a = n j + 1`, the same condition reads

```
    ⌊ d · a / e ⌋  mod  n   ∈   {1, …, n−2},    or  = n−1 with e | d·a.
```

So `⌊d·a/e⌋` is the Beatty sequence of slope `d/e`, and an opening index is a point at which the
Beatty sequences of the other disparities simultaneously avoid two residues modulo `n`. The flank
itself is free, since `d = e` gives `⌊a⌋ = a ≡ 1 (mod n)`. The clause `= n−1 with e | d·a` is the
closed endpoint of the band and cannot be dropped.

A certificate for whole-number disparities is checked in `n − 1` modular reductions and carries no
trace of how it was found. Commensurable real disparities are whole numbers after a common scaling,
which changes no loneliness, so Theorem 10 decides every commensurable subject.

**Every opening is a moment of global isolation.** The opening `(j + 1/n)/e` equals `a/(n e)` with
`a = n j + 1`, and `gcd(a, n) = 1`. By Theorem 4 this is a moment of global isolation of the arithmetic
progression whose common difference is `e`: at the gap between the subject and its flank, the
progression of that gap would isolate every runner at once. A certificate asks that the actual
runners, read against that progression, all stay clear. The opening index is always a unit of `ℤ/n`
— the class `a ≡ 1` — and no witness outside the units ever occurs.

---

## 7. The cascade: constructing the certificate

Theorems 9 and 10 say what loneliness *is*. The cascade produces a certificate, deterministically,
from the bound.

> **Definition.** A runner at disparity `d` is **in the way** at `t` when `‖d·t‖ < 1/n`. Its **next
> opening** after `t` is the least `(N + 1/n)/d` at or after `t`. The **move** sends `t` to the
> largest next opening among the runners in the way at `t`, and fixes `t` if none is. The
> **cascade** is the sequence of moves from the bound `1/(n·d_min)`.

A runner in the way at `t` stays in the way until its own next opening; the next opening is strictly
later than `t`; and each move lands on the opening of one runner in the way. There is no candidate
list: the runners in the way at `t` fix the move.

It is a recursion on the scale. At `1/(nD)` a runner of reading `s = d/D` is in the way exactly when
`|s − nN| < 1` for some `N`, and its next opening is `1/(n·d/(nN+1))`. So

```
    D₀ = d_min,    D_{j+1} = D_j · min{ s/(nN+1) : |s − nN| < 1 },    instant = 1/(n·D_final),
```

stopping when no reading lies within 1 of a multiple of `n`. Every factor is below 1, so the scales
decrease.

### 7.1 A complete example

Take the arrangement `1, 2, 3, 4, 5, 6, 8, 9`, `n = 8`, and the subject at speed `1`, whose
disparities are

```
    ds = { 1, 2, 3, 4, 5, 7, 8 },      d_min = 1,      bound t₀ = 1/8.
```

**Step 0, `t = 1/8`.**

| `d` | 1 | 2 | 3 | 4 | 5 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| `d·t` | 1/8 | 1/4 | 3/8 | 1/2 | 5/8 | 7/8 | 1 |
| `‖d·t‖` | 1/8 | 1/4 | 3/8 | 1/2 | 3/8 | 1/8 | **0** |

Runner 8 has completed exactly one relative lap and sits on the subject; it is the only runner in
the way. Its next opening is `(1 + 1/8)/8 = 9/64`. Runners 1 and 7 sit at distance exactly `1/8`, on
the boundary, and are clear — had the inequality been strict, this instant would have three
obstructions rather than one.

**Step 1, `t = 9/64`.**

| `d` | 1 | 2 | 3 | 4 | 5 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| `‖d·t‖` | 9/64 | 9/32 | 27/64 | 7/16 | 19/64 | **1/64** | 1/8 |

Runner 8 is now on the boundary. Runner 7, which sat on the boundary at step 0, reads `63/64` and is
in the way. Its next opening is `(1 + 1/8)/7 = 9/56`.

**Step 2, `t = 9/56`.**

| `d` | 1 | 2 | 3 | 4 | 5 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| `‖d·t‖` | 9/56 | 9/28 | 27/56 | 5/14 | 11/56 | **1/8** | 2/7 |

Every distance is at least `1/8`, the move fixes `t`, and the subject is lonely at `t = 9/56`.

**Reading off the certificate.** `9/56 = (1 + 1/8)/7` is the opening `j = 1` of the flank `e = 7`, so
the certificate is `(7, 1)`, checked in seven modular reductions.

**The rules of §10 reach this set by another route.** Its largest disparity exceeds seven times its
least, so R1 does not apply; one disparity is divisible by 8, so R3 does not. No disparity is
divisible by 6, so R6 applies at `q = 6` and names the instant `t = 1/6`, at which the distances are
`1/6, 1/3, 1/2, 1/3, 1/6, 1/6, 1/3`, all at least `1/8`. That instant is not an opening; Theorem 9
moves it back to one. The cascade always terminates when the subject is lonely and needs no
condition; a rule bypasses the construction but applies only where its condition holds.

### 7.2 Construction and Localisation

> **Theorem 11 (Construction).** If some step of the cascade is fixed, at `T`, then: (a) the subject
> is lonely at `T`; (b) no positive instant before `T` is lonely; (c) `T` is the end of a run of
> lap-completion zones from the bound, each containing the end of the one before, and every such run
> that ends at a lonely instant ends at `T`; and (d) some runner sits at exactly `1/n`.

*Proof.* **(a)** "No runner in the way at `t`" says `‖d·t‖ ≥ 1/n` for every disparity.

**(b)** Nothing below the bound is lonely, by Theorem 7. Let `t` be a step that is not fixed, with
image `T'`. Every runner in the way at `t` remains in the way until its own next opening, and `T'` is
the **largest** of those, so every instant in `[t, T')` has some runner in the way. Chaining along the
finitely many steps from the bound to `T` leaves no lonely instant below `T`.

**(c)** Each move lands on the far endpoint of one obstructing runner's zone
`((N − 1/n)/d, (N + 1/n)/d)`, and the runner that sets the following move is in the way there, so
each zone contains the end of the one before. For the second clause, the open zones of such a run
cover every instant from the bound up to its end. If it ended after `T`, then `T` would lie inside
one of its zones and not be lonely; if it ended before `T`, its end would be a lonely instant before
`T`, against (b).

**(d)** If `T` is the bound then `d_min` reads exactly `1/n` there. Otherwise `T` is the next opening
`(N + 1/n)/d` of some runner `d`, which then sits at exactly `1/n`. ∎

So a fixed point is a lonely instant, and it is always an opening — which is why the certificate has
the form `(e, j)`. Part (b) says it is also the first; nothing in Layer 2 uses that.

> **Theorem 12 (Cascade localisation).** For every positive instant `s` at which the subject is
> lonely: `s ≥ 1/(n·d_min)`; every step of the cascade is at most `s`; `s` is matched by a lonely
> opening no later than it; and **the cascade reaches a fixed point**.

*Proof.* The first clause is Theorem 7; the third is Theorem 9.

For the second, induct along the cascade. Suppose a step `t ≤ s` is not fixed, and let `d` be in the
way at `t`, so `d·t ∈ (N − 1/n, N + 1/n)`. Since `d·s ≥ d·t` and `d` is clear at `s`, `d·s ≥ N + 1/n`,
so `d`'s next opening is at most `s`. The move goes to the largest of those, so its image is at most
`s`.

For the last, every step lies in the finite set of openings `(k + 1/n)/d` with `0 ≤ k ≤ ⌈d·s⌉`, and a
sequence that strictly increases until it is fixed cannot remain in a finite set without being fixed.
∎

> **Theorem 13.** For every injective arrangement `v` and every subject `k`, with no hypothesis
> beyond injectivity, the following are equivalent:
>
> 1. the subject is lonely at some positive instant;
> 2. the cascade reaches a fixed point;
> 3. the cascade's instants are bounded above;
> 4. some flank has a lonely opening;
> 5. each runner is in the way at only finitely many steps.

*Proof.* **(2) ⟹ (1)** is Theorem 11(a), and **(1) ⟹ (2)** is Theorem 12.
**(2) ⟹ (3)**: a cascade that reaches a fixed point is constant from then on.
**(3) ⟹ (2)**: bounded steps lie in a finite set of openings, and a strictly increasing sequence
cannot stay in a finite set.
**(1) ⟺ (4)** is Theorem 9.
**(2) ⟹ (5)**: beyond a fixed point no runner is in the way.
**(5) ⟹ (2)**: there are finitely many runners, so beyond some step none is in the way. ∎

Nothing in the argument mentions the order except through the threshold `1/n`.

---

## 8. The framework stated

Collecting Theorems 9–13:

```
    loneliness   ⟺   lonely opening   ⟺   certificate   ⟺   cascade halts
```

Each arrow is an equivalence, proved in both directions, with no hypothesis beyond injectivity of the
speeds, and uniformly in the order. For a fixed subject of a fixed arrangement:

* the subject is lonely at some instant **exactly when** one of its own disparities has an opening at
  which every other disparity is clear (Theorem 9);
* for whole-number disparities that condition is **exactly** a finite conjunction of modular
  memberships, decided by the `e` indices `0 ≤ j < e` (Theorem 10);
* the cascade, started at the bound and driven only by the runners in the way, halts **exactly when**
  the subject is lonely, with a certificate readable off the halting state
  (Theorems 11–13).

**This fixes what a proof is.** Because the equivalence is unconditional and runs both ways,
exhibiting a certificate is not evidence for the conjecture at that configuration; it **is** the
conjecture there, written in the framework's own objects. A case discharged is proved outright, with
nothing further owed on it, and a proof of the whole is a finite list of such cases whose conditions
are jointly exhaustive.

Three shapes ask for more than that, and none is used below.

| shape | why it is more than is owed |
|---|---|
| a **cap** — a bound on how far the cascade runs | the conjecture asks for loneliness at *some* instant; bounding how far a run goes is an obligation it does not carry |
| a **move count** — "the `N`-th step is fixed" | a bound on *how many* steps, where the conjecture places none |
| a **uniform theorem** — one statement or formula covering every disparity set | the form is conditions with proofs; a closed form for `e` or `j` is more than a case requires |

---

# LAYER 2 — CONDITIONS FOR A CERTIFICATE

---

## 9. The obligation at order `n`

Layer 1 settles, for a fixed configuration, what loneliness is:

```
    ∃ t > 0 lonely      ⟺      ∃ (e, j) with the opening (j + 1/n)/e lonely.
```

The conjecture at order `n` is the statement that such an `(e, j)` exists for every subject of every
arrangement. For whole-number disparities it can be put on sets.

> **Lemma 14 (the obligation on sets).** The conjecture at order `n` holds for integer speeds if and
> only if every set of at most `n − 1` positive integers is lonely at threshold `1/n`.

*Proof.* The disparities of a subject form such a set. Conversely, given a set `S` of at most `n − 1`
positive integers, choose further distinct positive integers so that `S` and them make `n − 1` in
all, and take the arrangement consisting of `0` and those speeds. The subject at `0` is lonely at
some `t` by hypothesis, and the disparities in `S` are among its own, so `S` is lonely at `t`. ∎

So a smaller set is never harder than a larger one at the same order, and an induction on the number
of disparities at a fixed order is available. Everything in Layer 2 works inside this obligation.

---

## 10. Certificate rules

A rule is a condition on a disparity set together with a proof that a certificate exists whenever the
condition holds. Where a rule's condition holds, the conjecture is established for that
configuration. All seven hold at every order.

| rule | condition on the disparity set | what it exhibits |
|---|---|---|
| **R1** the ratio condition | `d_max ≤ (n−1)·d_min` | the bound `1/(n·d_min)`, the opening `j = 0` of `d_min` |
| **R2** the band ladder | every disparity in `[(nL+1)e/(nj+1), (nL+n−1)e/(nj+1)]` for some `j, L` | the `j`-th opening of `e` |
| **R3** the non-resonant flank | `gcd(e, n) = 1` and no disparity divisible by `n` | the opening `a = e·s`, `e·s ≡ 1 (mod n)` |
| **R4** the divisible flank | `n \| e`, and either (a) `e` largest, no other disparity divisible by `n`, fewer than `φ(n)` of the others coprime to `n`; or (b) every other multiple of `n` clear at `j = 0`, and the others of total cost `∑ max(2, gcd(d, n)) < n` | the opening `a = k·e + 1`, `k < n` |
| **R5** the descent | the disparities not divisible by `n` have total cost `∑ max(2, gcd(e, n)) < n` | reduces to the cofactors of the multiples, same order |
| **R6** the divisibility rules | for some `2 ≤ q ≤ n`, no disparity divisible by `q` | the instant `t = 1/q` |
| **R7** the free disparity | one disparity off the `P`-lattice of the rest; the rest at most three and sharing it | a lonely instant, `n ≥ 4` |

R1–R4 exhibit an opening directly. R5, R6 and R7 exhibit a lonely instant, which Theorem 9 converts
to a certificate. R1–R6 are stated for whole-number disparities; R7 for real ones.

### 10.1 Scale control — R1

This is Theorem 8.

### 10.2 Band control — R2

At the opening `j` of `e`, the reading of `d` is `(d/e)(nj+1)`, clear when it lies between
consecutive multiples of `n` at distance at least 1 from each. Requiring **one shared lap counter `L`
for every disparity** gives `(nL+1)e/(nj+1) ≤ d ≤ (nL+n−1)e/(nj+1)`: all disparities in a single
band. At `j = L = 0` with the flank at `d_min` the band is `[d_min, (n−1)d_min]`, which is R1. Eliminating
`L` gives a necessary condition, `(d_max − d_min)(nj+1) ≤ e(n−2)`.

### 10.3 Non-resonant structure — R3

Let `gcd(e, n) = 1` and no disparity be divisible by `n`. Take `s` with `e·s ≡ 1 (mod n)` and
`a = e·s`, which is `≡ 1 (mod n)` and so an opening index. The reading `d·a mod n e` is
`e·(d·s mod n)`, in the band `[e, (n−1)e]` exactly when `n ∤ d·s`, and since `s` is a unit that is
exactly `n ∤ d`. No primality is used: the hypothesis is that the flank is **non-resonant**. R3 is
the `q = n` instance of R6 with an opening named.

### 10.4 Resonant structure — R4 and R5

**R4, the divisible flank.** Let `n | e`, let `e` exceed every other disparity, let no other
disparity be divisible by `n`, and let fewer than `φ(n)` of them be coprime to `n`.

*Proof.* For `k` a unit of `ℤ/n` take `a = k·e + 1`; since `n | e`, `a ≡ 1 (mod n)`. The flank reads
`e`, in the band. Another disparity `d < e` reads `d·a mod n e = e·(d·k mod n) + d`, which lies in
`[e, (n−1)e]` exactly when `d·k mod n ∉ {0, n−1}`. Since `k` is a unit and `n ∤ d`, `d·k ≢ 0`. If
`d·k ≡ −1 (mod n)` then `n | d·k + 1`, so `gcd(d, n) = 1`: a disparity sharing a factor with `n` rules
out no unit. A disparity coprime to `n` rules out exactly one unit, since `k ↦ d·k mod n` is
injective. There are `φ(n)` units and fewer than `φ(n)` coprime disparities, so some unit survives.
∎

At a prime order every non-multiple is coprime and `φ(n) = n − 1`, so for a subject with its full
`n − 2` non-flank disparities the budget condition holds automatically. At `n = 8` it asks that at
most three of the non-flank disparities be odd.

**R4 without the largest hypothesis.** Let `n | e`, with no condition on the size of `e`. Suppose
every other multiple of `n` is clear at the opening `j = 0` of `e`, and the disparities not divisible
by `n` have total cost `∑ max(2, gcd(d, n))` below `n`.

*Proof.* Take `a = k·e + 1` with `0 ≤ k < n`, so that `a ≡ 1 (mod n)` and `a < n·e`. The flank reads
`e`. A multiple `d = n·m` reads `d + (n·e)(m·k) ≡ d (mod n·e)`, the same as at `j = 0`, so it is
clear by hypothesis. A non-multiple `d = q·e + r`, with `q = ⌊d/e⌋`, has `r ≥ 1` (otherwise `e | d`,
and so `n | d`), and reads `e·((d·k + q) mod n) + r`. That lies in `[e, (n−1)e]` exactly when
`d·k + q ≢ 0, −1 (mod n)`. Each congruence `d·k ≡ c (mod n)` has at most `gcd(d, n)` solutions below
`n`, since two solutions differ by a multiple of `n/gcd(d, n)`; and when `gcd(d, n) ≥ 2` the residues
`−q` and `−q − 1` cannot both be solvable, since `gcd(d, n)` would divide both and hence `1`. So `d`
forbids at most `max(2, gcd(d, n))` of the `n` values of `k`, the costs add up to less than `n`, and
some `k` is forbidden by none. ∎

When `e` is the largest disparity, `q = 0` for every other disparity and form (a) is usually the
sharper, since restricting `k` to units removes the residue `0` for free. Form (b) needs no ordering
at all.

**R5, the descent.** Split the disparities into a set `E` of disparities not divisible by `n` and
the multiples `n·m`, `m ∈ M`. Give each `e ∈ E` the **cost** `max(2, gcd(e, n))`.

> **Theorem 15 (Descent).** If the costs of `E` add up to less than `n`, then `E ∪ n·M` is lonely at
> threshold `1/n` if and only if `M` is lonely at the same threshold `1/n`.

*Proof.* If `t` is lonely for `E ∪ n·M`, then `τ = n·t` clears every `m ∈ M`, since `m·τ = (n·m)·t`.

Conversely let `τ > 0` clear `M`, and try `t = (τ + k)/n` for `k = 0, 1, …, n − 1`. Each multiple
reads `(n·m)·t = m·τ + m·k`, unchanged mod 1, so it stays clear. A disparity `e ∈ E` reads
`x + e·k/n` with `x = e·τ/n`; call the shift `k` *blocked* for `e` when that reading is within `1/n`
of an integer `K_k`. For two blocked shifts `k, k'`, the integer `σ(k, k') = e(k − k') − n(K_k − K_{k'})`
is `n` times a difference of two quantities each below `1/n` in absolute value, so `|σ(k, k')| ≤ 1`;
and the steps add, `σ(k, k'') = σ(k, k') + σ(k', k'')`.

* If `gcd(e, n) = 1`, a step `σ(k, k') = 0` gives `n | e(k − k')`, hence `k = k'`. Three distinct
  blocked shifts would give three nonzero steps of size 1, one the sum of the other two, which is
  impossible. So at most two shifts are blocked.
* If `g = gcd(e, n) ≥ 2`, then `g` divides every step, so every step is `0`, and `n | e(k − k')` gives
  `(n/g) | (k − k')`. All blocked shifts lie in one residue class mod `n/g` below `n`, so there are at
  most `g`.

So `e` blocks at most `max(2, gcd(e, n))` of the `n` shifts. The costs add up to less than `n`, so some
shift is blocked for no `e ∈ E`, and at that `t` every disparity is clear. ∎

A multiple of `n` would cost `n` on its own, so the condition already excludes one from `E`. With `E`
empty the theorem says that dividing every disparity by `n` changes no loneliness. With `E = {e}` and
`n ≥ 3` the condition always holds, so a single non-multiple can always be removed. At `n = 8` a
disparity costs 4 if it is `≡ 4 (mod 8)` and 2 otherwise, so up to three non-multiples not `≡ 4`, or
one `≡ 4` with one other, can be removed at once. R5 therefore replaces a set by one with fewer
disparities, or with the same number and smaller ones, at the same order. By Lemma 14 that is the
same obligation, not a lower-order conjecture.

### 10.5 Divisibility — R6

For `2 ≤ q ≤ n`, the instant `1/q` is lonely exactly when no disparity is divisible by `q`: the
reading of `d` is `(d mod q)/q`, at distance at least `1/n` from the integers precisely when it is
nonzero, since `q ≤ n`. That is `n − 1` rules, each an equivalence. Their instants depend on the
disparity set only through divisibility and are indexed by `q`, not by any member of the set.

### 10.6 The free disparity, and the real case — R7

Suppose one disparity `x` sits off the lattice `{ y : P·(y/f) ∈ ℤ }` that the others share for some
integer `P` — in particular if `x` has irrational ratio to them — and the others number at most three.
Then, at `n ≥ 4`, the subject is lonely: the others are settled at some flank, and `x`, being off the
lattice, cannot obstruct a whole residue class of that flank's openings, so one of them is clear of
`x` too. A disparity incommensurable with the others therefore helps. The sets that resist the rules
are the commensurable ones.

### 10.7 The units of `ℤ/n`

The units of `ℤ/n` appear throughout, and primality nowhere. They index the moments of global
isolation of the progression (Theorem 4); every opening index `a = n j + 1` is one (§6); R3's multiplier `s` and R4's
multiplier `k` are units, and R4's budget is their number, `φ(n)`. A general witness has the same
shape: any `a ≡ 1 (mod n)` with `a < n·e` at which every reading `d·a mod n·e` lies in `[e, (n−1)e]`
is a certificate.

### 10.8 The failure conditions

A rule whose failure is characterised is a case in a case analysis, and cases compose.

> **R6 fails** exactly when every modulus in `[2, n]` divides some disparity.
>
> **R3 applies** exactly when no disparity is divisible by `n` and some disparity is coprime to it.
>
> **R5 applies** exactly when the disparities not divisible by `n` have total cost below `n`; at
> `n ≥ 3` this includes every set with at most one such disparity.

R3's condition is the `q = n` instance of R6, which does not need a coprime flank to exist.

---

## 11. How the rules compose

Write `M` for the number of disparities divisible by `n`, and call the total cost of the others,
`∑ max(2, gcd(e, n))`, their **cost**.

| region | standing |
|---|---|
| `M = 0` | **settled** — R6 at `q = n` names `t = 1/n` |
| cost below `n` | **passed down** — R5 reduces the set to its cofactors, at the same order |
| `M ≥ 1` and cost at least `n` | the middle |

> **Corollary 16 (composition).** Let `n ≥ 3`. If the disparities of a set not divisible by `n` have
> cost below `n`, and the set of cofactors of its multiples of `n` satisfies any condition on the
> list — or is itself passed down by R5 to a set that does — then the set is lonely at threshold
> `1/n`.

*Proof.* By Theorem 15, applied as many times as the passing down continues. Each application leaves
fewer disparities, or the same number with smaller ones, so it stops. ∎

So R5 is not an eighth region but a way of carrying every other rule into sets it would not reach
directly. For example, at `n = 5` the set `{1} ∪ 5·{1, 2, 3} = {1, 5, 10, 15}` is passed down to
`{1, 2, 3}`, which R1 settles.

The middle is not the region that remains. R1, R2, R4, R6 at `q < n` and R7 are conditions of other
kinds and reach into it without reference to the cost. Both middle examples in this paper are settled:
§7.1's set by R6 at `q = 6`, and §12's `{1, 2, 3, 5}` by R4. What remains is the part of the middle on
which every rule fails.

---

## 12. The middle

A disparity divisible by `n`, say `n·m'`, reads

```
    (n·m'·a)  mod  n e   =   n · ( (m'·a)  mod  e ),
```

so its band condition is a band condition at modulus `e`, and it depends on `a` only through
`ρ = a mod e`. A disparity not divisible by `n` needs all of `a`. So the middle is a copy of the same
problem one scale down, constrained additionally by the non-multiples:

```
    middle   =   resonant component at modulus e   +   non-resonant constraints.
```

R3 is the case where the first part is empty. R5 removes the non-multiples whenever their cost is
below `n`: the shift by `k/n` moves only the non-multiples, each blocks at most `max(2, gcd(e, n))` of
the `n` shifts, and a free shift remains. When the cost reaches `n`, that count no longer guarantees a
free shift, and a further condition is needed.

**Example.** Take `{1, 2, 3, 5}` at order five. The resonant part is `{5}`; `1, 2, 3` are
non-multiples. The set is lonely, and R4 reaches it: `5` is the largest disparity and the only
multiple of `5`, and the three others are coprime to `5` with `3 < φ(5) = 4`. It names the flank
`e = 5` at `j = 1`, `a = 6`: the readings `6, 12, 18, 5` modulo 25 all lie in `[5, 20]`. The flank
`e = 3` at `j = 2` also clears.

**What the middle needs.** Conditions on the disparity set under which some flank and some index can
be shown to exist — by pigeonhole, induction, contradiction or construction.

---

## 13. Further structure

The development record, Palelei (2026b), contains further machine-checked results on the certificate
condition: its reparameterisation by the cascade's contraction, the candidate family and its
rescaling, a reflection symmetry and its two obstructing residues, an edge-witness property of
finite systems of congruence-interval constraints, and the geometry of the clear window. None is
needed for the results above.

---

## 14. What is proved, and what remains

**Proved — the characterisation (§§2–8).** Loneliness is a two-neighbour condition; global
isolation is the regular `n`-gon, and a moment of global isolation survives a change of speeds exactly
under Theorem 6; the number of simultaneously lonely runners is never `n − 1`. Every lonely instant is
represented by a lonely opening, a moment of global isolation at the runner's own scale, and
conversely; for whole-number disparities the opening condition is a finite conjunction of modular
memberships; the cascade halts exactly when the subject is lonely. All of this holds at every order
with no hypothesis beyond injectivity.

**Proved — the cases (§§10–11).** Seven rules, each a condition with a proof that a certificate
exists under it, at every order. Every configuration satisfying any of them is settled, and at every
order there are infinitely many. R5 carries every rule into further sets by passing them down to
their cofactors (Corollary 16).

**Remaining.** Conditions covering the disparity sets in the middle on which every rule fails, each
discharged by showing a certificate exists. When they are found, the isolation of the ideal is shown
to be shifted and never destroyed at every order. The list is incomplete; the entries on it are not
provisional.

---

## 15. General `n`

Nothing in §§2–8 depends on the order except through the threshold `1/n`. The rules of §10 and the
trichotomy of §11 hold at every order, with the stated conditions (R5's cost condition, which forces
`n ≥ 3` when a non-multiple is present; `n ≥ 4` for R7). Two further facts are proved at general `n`:

* subjects with at most three distinct disparities are lonely at every order `n ≥ 4`, and subjects
  with at most four at every order `n ≥ 6`;
* the fastest and the slowest runner of every arrangement carry the full `n − 1` distinct
  disparities, since every difference has one sign there, so the full-cardinality case cannot be
  avoided.

---

## 16. Verification

The results below are formalised in Lean 4 with mathlib, with no `sorry`, no `native_decide`, no
`unsafe` and no `admit`, and every theorem depends on at most `propext`, `Classical.choice` and
`Quot.sound`. The formalisation verifies the stated claims; it is not part of the mathematical
argument.

| result | declaration |
|---|---|
| Theorem 1 | `lonely_iff_flanks` |
| Theorem 2 | `gaps_all_eq` |
| Corollary 3 | `lonely_ne_pred` |
| Theorem 4 *(Palelei 2026a)* | `ap_reference` (the instants are global), `card_nonresonant_eq_totient` (there are `φ(n)` classes); the converse is Palelei (2026a) |
| Theorem 5 | `displacement` |
| Theorem 6 | `survives_of` (forward), `int_of_pairwise_sep` and `not_dvd_of_sep` (converse), `commensurable` |
| Theorem 7 | `bound_le_of_lonely` |
| Theorem 8 | `ratio_condition` |
| Theorem 9 | `lonely_iff_exists_flank` |
| Theorem 10 | `clear_at_int_flank_iff_emod`, `flankWorks_iff`, `band_iff_beatty` |
| §6, openings are units | `witness_is_unit` |
| Theorem 11 | `lonely_of_fixed`, `first_isolation_of_fixed`, `first_isolation_as_zone_run`, `run_end_eq_fixed_point`, `flank_at_first_isolation` |
| Theorem 12 | `bound_le_of_lonely`, `iterate_le_of_lonely`, `exists_fixed_of_lonely` |
| Theorem 13 | `doors_agree` (general `n`) |
| R1 | `ratio_condition`, `clear_at_last_opening` |
| R2 | `lonely_flank_of_band`, `bandFits_forces_spread`, `spread_zero_iff_ratio` |
| R3 | `flankWorks_of_coprime`, `flankWorks_of_not_resonant` |
| R4 | (a) `flankWorks_of_dvd_totient`, `exists_goodK_totient`; (b) `flankWorks_of_dvd_budget`, `card_forbidden_le`, `card_congr_le` |
| R5, Theorem 15 | `lonely_multi_descent_iff`, `card_bad_le`, `bad_step`; one non-multiple: `lonely_descent_iff`, `exists_clear_shift`; all multiples: `lonely_scale_iff` |
| R6 | `ratWorks_le_order`, `ratWorks_one_iff`, `ratWorks_family_fails_iff` |
| R7 | `lonely_of_one_free`, `lonely_of_one_irrational` |
| the general witness form | `flankWorks_of_witness` |
| §10.8 | `multiples_trichotomy`, `self_rule_applies_iff`, `coprime_rule_applies_iff`, `divisibility_subsumes_coprime_flank` |
| §15 | `exists_lonely_of_card_le_three`, `exists_lonely_of_card_le_four`, `diffSet_card_eq_of_extreme` |

Lemma 14 and Corollary 16 are proved in the text from the declarations above.

---

## 17. Glossary

**arrangement** — an injective assignment of real speeds to the `n` runners.
**band** — the readings at distance at least 1 from `n·ℤ`; equivalently, circle distance at least
`1/n`.
**boundary set method** — evaluate open cells and boundary points separately; never merge touching
zones; decide a cell by an interior point.
**cascade** — the sequence of moves from the bound.
**certificate** — a flank and an opening index `(e, j)` whose opening is lonely.
**disparity** — one of the subject's speed differences `|vᵢ − v_k|`.
**flank** — a disparity used as the denominator of an opening.
**cost** — of a disparity `e` not divisible by `n`, the number `max(2, gcd(e, n))`.
**global isolation** — all `n` runners lonely at the same moment; for equally spaced speeds it occurs
`φ(n)` times per period.
**in the way** — at circle distance strictly less than `1/n` from the subject.
**move** — the step sending `t` to the largest next opening among the runners in the way.
**non-resonant** — coprime to `n`.
**opening** — the instant `(j + 1/n)/e`, at which the flank `e` sits at exactly `1/n`.
**reading** — the quantity `n·d·t`, or `s·u` in the local coordinate.
**resonant** — sharing a factor with `n`.
**subject** — the runner whose loneliness is at issue.

---

## Data and code availability

The Lean 4 development accompanies this paper as supplementary material. The development record,
declaration indices and further results are in Palelei (2026b).

---

## References

Barajas, J. and Serra, O. (2008). The lonely runner with seven runners. *Electron. J. Combin.* 15,
R48.

Bohman, T., Holzman, R. and Kleitman, D. (2001). Six lonely runners. *Electron. J. Combin.* 8(2),
R3.

Palelei, T. (2026a). The φ(n) law for arithmetic progressions in the lonely runner conjecture:
algebraic structure and computational fragility. Zenodo. https://doi.org/10.5281/zenodo.18158886

Palelei, T. (2026b). A bounding framework for eight lonely runners: development record, statements
and verification apparatus. Technical report, supplementary material to the present paper.

Rosenfeld, M. (2025). The lonely runner conjecture holds for eight runners. arXiv:2509.14111.

Sungkawichai, T. and Trakulthongchai, T. (2026). Eleven, twelve, and thirteen lonely runners.
arXiv:2604.23906.
