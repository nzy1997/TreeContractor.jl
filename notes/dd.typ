#import "@preview/clear-iclr:0.7.0": *
#import "@preview/cetz:0.4.2": canvas, draw

#show: iclr.with(
  title: [Rank-Width for \#SAT: A Tensor Network Perspective],
  authors: (
    (
      names: ([Anonymous],),
      affilation: [Tutorial Notes],
      address: [],
      email: "cacate0129@gmail.com",
    ),
  ),
  keywords: ("rank-width", "SAT", "tensor networks", "treewidth"),
  abstract: [
    We provide a pedagogical introduction to rank-width and its application to the propositional model counting problem \#SAT. We bridge the gap between rank-width and familiar concepts from tensor network theory and treewidth-based algorithms. Through a concrete example, we demonstrate how rank-width can lead to exponential speedups over treewidth-based approaches for certain graph families.
  ],
  accepted: none,
  bibliography: none,
)

= Introduction

The paper by Ganian, Hliněný, and Obdržálek @ganian2010 provides a parameterized polynomial algorithm for \#SAT with runtime single-exponential in the _rank-width_ of a formula. This improves upon previous algorithms based on _clique-width_, since rank-width can be exponentially smaller.

For practitioners familiar with tensor networks and treewidth, rank-width offers an alternative decomposition strategy that captures a fundamentally different notion of "width" --- measuring the *rank of certain matrices* rather than the size of separators.

= Background: Treewidth vs Rank-Width

#figure(
  table(
    columns: 3,
    stroke: none,
    table.hline(),
    table.header([*Concept*], [*Treewidth*], [*Rank-Width*]),
    table.hline(),
    [Decomposition], [Tree decomposition (bags)], [Rank decomposition (binary tree)],
    [Width measure], [max bag size $- 1$], [max rank of cut matrices],
    [Algorithm complexity], [$O(2^(t w) dot n)$ for \#SAT], [$O(2^(r w) dot n)$ for \#SAT],
    [Relationship], [$r w(G) <= t w(G) + 1$], [Can be exponentially smaller],
    table.hline(),
  ),
  caption: [Comparison of treewidth and rank-width approaches.],
) <tab:comparison>

The key relationship is:
$ r w(G) <= t w(G) + 1 $
but rank-width can be *exponentially smaller* than treewidth for certain graph families.

= Rank Decomposition

== Definition

A *rank decomposition* of a graph $G = (V, E)$ consists of:
1. A *binary tree* $T$ with leaves bijecting to vertices $V$
2. For each edge $e$ of $T$, removing $e$ partitions the leaves into two sets $A$ and $overline(A) = V without A$

== The Cut Matrix

For a partition $(A, overline(A))$, define the *cut matrix* $M_A$ over $bb(F)_2$ (the binary field):
- Rows indexed by $A$, columns by $overline(A)$
- $M_A [u, v] = 1$ if $(u, v) in E$, else $0$

The *rank-width* is:
$ r w(G) = min_T max_(e in T) op("rank")_(bb(F)_2) (M_A) $

// == Tensor Network Interpretation

// If we think of the adjacency matrix as a tensor, the cut matrix is exactly what we get when we *partition the indices* and reshape into a matrix. The rank over $bb(F)_2$ measures how "entangled" the two sides are --- analogous to the bond dimension in tensor network terminology.

= A Concrete \#SAT Example: Dense Formulas

To see the power of rank-width, we need a formula where rank-width is *much smaller* than treewidth. This happens with *dense* formulas where every variable appears in every clause.

== A Dense SAT Formula

Consider a formula where *every variable appears in every clause*:
$ psi = and.big_(i=1)^4 C_i $
where each $C_i$ contains all 4 variables $x_1, x_2, x_3, x_4$ (with various polarities).

For concreteness:
- $C_1 = (x_1 or x_2 or x_3 or x_4)$
- $C_2 = (not x_1 or x_2 or x_3 or x_4)$  
- $C_3 = (x_1 or not x_2 or x_3 or x_4)$
- $C_4 = (x_1 or x_2 or not x_3 or x_4)$

== The Incidence Graph is $K_(4,4)$

Since every variable appears in every clause, the incidence graph is the *complete bipartite graph* $K_(4,4)$:

#figure(
  canvas(length: 1cm, {
    import draw: *
    
    let var-color = rgb("#4a90d9")
    let clause-color = rgb("#e07b53")
    
    // Variables on left
    for (i, label) in ((0, $x_1$), (1, $x_2$), (2, $x_3$), (3, $x_4$)) {
      circle((0, -i * 0.9), radius: 0.3, fill: var-color.lighten(70%), stroke: var-color, name: "x" + str(i+1))
      content((0, -i * 0.9), text(8pt, label))
    }
    
    // Clauses on right
    for (i, label) in ((0, $C_1$), (1, $C_2$), (2, $C_3$), (3, $C_4$)) {
      rect((2.7, -i * 0.9 - 0.2), (3.3, -i * 0.9 + 0.2), fill: clause-color.lighten(70%), stroke: clause-color, name: "C" + str(i+1))
      content((3, -i * 0.9), text(8pt, label))
    }
    
    // All edges (complete bipartite)
    set-style(stroke: gray.lighten(30%) + 0.5pt)
    for i in range(4) {
      for j in range(4) {
        line("x" + str(i+1), "C" + str(j+1))
      }
    }
  }),
  caption: [The incidence graph $K_(4,4)$: every variable connects to every clause.],
) <fig:k44>

== Treewidth vs Rank-Width of $K_(n,n)$

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([Graph], [Treewidth], [Rank-width]),
    [$K_(4,4)$], [$4$], [$1$],
    [$K_(n,n)$], [$n$], [$1$],
    [$K_(100,100)$], [$100$], [$1$],
  ),
  caption: [Complete bipartite graphs have constant rank-width but linear treewidth.],
) <tab:knn>

*Why is the rank-width only 1?*

The biadjacency matrix of $K_(n,n)$ is the all-ones matrix $J$:
$ M = mat(1, 1, 1, 1; 1, 1, 1, 1; 1, 1, 1, 1; 1, 1, 1, 1) $

Over $bb(F)_2$, this has *rank 1* because all rows are identical! We can write:
$ M = vec(1, 1, 1, 1) dot mat(1, 1, 1, 1) = bold(1) dot bold(1)^top $

== The Dramatic Speedup

For the dense formula $psi$ with $n$ variables and $n$ clauses:
#v(20pt)

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([Approach], [State space size], [For $n = 100$]),
    [Treewidth-based], [$2^n$], [$2^(100) approx 10^(30)$],
    [Rank-width-based], [$2^1 = 2$], [$2$],
  ),
  caption: [Exponential difference in state space for dense formulas.],
) <tab:speedup>

The rank-width algorithm only needs to track *2 equivalence classes* instead of $2^(100)$ configurations!

== Why Does This Work? (SVD Perspective)

For tensor network practitioners, this is best understood via *SVD/low-rank factorization*.

The biadjacency matrix $M = bold(1) bold(1)^top$ factors as:
$ M = underbrace(vec(1, 1, 1, 1), U) dot underbrace(mat(1, 1, 1, 1), V^top) $

For any cut partitioning variables into $A$ and $overline(A)$, a partial assignment $bold(v) in bb(F)_2^(|A|)$ produces a boundary effect:
$ bold(v)^top M = (bold(v)^top bold(1)) dot bold(1)^top $

The term $(bold(v)^top bold(1)) in bb(F)_2$ is just the *parity* (sum mod 2) of the assignment bits. This gives only 2 possible values:
- Parity 0 $arrow.r$ boundary effect $= bold(0)$
- Parity 1 $arrow.r$ boundary effect $= bold(1)$

#figure(
  canvas(length: 1cm, {
    import draw: *
    
    let box-fill = rgb("#f8f8f8")
    
    // Left: high-dimensional space
    rect((-0.5, -1), (1.5, 1), fill: box-fill, stroke: gray)
    content((0.5, 1.3), text(9pt)[$bb(F)_2^n$])
    content((0.5, 0), text(8pt)[$2^n$ partial\ assignments])
    
    // Arrow to middle
    line((1.7, 0), (2.8, 0), mark: (end: "stealth"), stroke: black + 0.8pt)
    content((2.25, 0.4), text(8pt)[$times bold(1)$])
    
    // Middle: low-dimensional space
    rect((3, -0.5), (4.5, 0.5), fill: rgb("#e8f4e8"), stroke: rgb("#4a4"))
    content((3.75, 0.8), text(9pt)[$bb(F)_2^1$])
    content((3.75, 0), text(8pt)[2 classes])
    
    // Arrow to right
    line((4.7, 0), (5.8, 0), mark: (end: "stealth"), stroke: black + 0.8pt)
    content((5.25, 0.4), text(8pt)[$times bold(1)^top$])
    
    // Right: boundary effect
    rect((6, -1), (8, 1), fill: box-fill, stroke: gray)
    content((7, 1.3), text(9pt)[$bb(F)_2^n$])
    content((7, 0), text(8pt)[boundary\ effects])
  }),
  caption: [SVD bottleneck: $2^n$ assignments compress to just 2 equivalence classes through the rank-1 factorization.],
) <fig:svd>

#v(20pt)
This is exactly the *bond dimension* in tensor networks: when you SVD a matrix of rank $r$, the bond dimension is $r$. Here $r = 1$, so we only need to track $2^1 = 2$ states instead of $2^n$.

== Connecting to Counting: Full DP Walkthrough

Let's trace through the algorithm step-by-step on a simpler instance to see actual numbers.

=== Setup: A 2×2 Dense Formula

Consider $phi = C_1 and C_2$ with variables $x_1, x_2$ where:
- $C_1 = (x_1 or x_2)$ — satisfied unless both false
- $C_2 = (not x_1 or not x_2)$ — satisfied unless both true

The incidence graph is $K_(2,2)$ (every variable in every clause).

#figure(
  canvas(length: 1cm, {
    import draw: *
    let var-color = rgb("#4a90d9")
    let clause-color = rgb("#e07b53")
    
    circle((0, 0), radius: 0.3, fill: var-color.lighten(70%), stroke: var-color, name: "x1")
    content((0, 0), text(8pt, $x_1$))
    circle((0, -1.2), radius: 0.3, fill: var-color.lighten(70%), stroke: var-color, name: "x2")
    content((0, -1.2), text(8pt, $x_2$))
    
    rect((1.7, -0.2), (2.3, 0.2), fill: clause-color.lighten(70%), stroke: clause-color, name: "C1")
    content((2, 0), text(8pt, $C_1$))
    rect((1.7, -1.4), (2.3, -1.0), fill: clause-color.lighten(70%), stroke: clause-color, name: "C2")
    content((2, -1.2), text(8pt, $C_2$))
    
    set-style(stroke: gray + 0.6pt)
    line("x1", "C1")
    line("x1", "C2")
    line("x2", "C1")
    line("x2", "C2")
  }),
  caption: [Incidence graph $K_(2,2)$ for $phi = (x_1 or x_2) and (not x_1 or not x_2)$.],
)

=== The Rank Decomposition Tree

#figure(
  canvas(length: 1cm, {
    import draw: *
    let node-fill = rgb("#f0f0f0")
    let leaf-fill = rgb("#e8f4e8")
    let cut-color = rgb("#e74c3c")
    
    // Root
    circle((0, 0), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "root")
    content((0, 0.5), text(8pt)[root])
    
    // Level 1
    circle((-1.5, -1.2), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "L")
    circle((1.5, -1.2), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "R")
    
    // Leaves
    circle((-2.2, -2.4), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "x1")
    content((-2.2, -2.4), text(8pt, $x_1$))
    circle((-0.8, -2.4), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "C1")
    content((-0.8, -2.4), text(8pt, $C_1$))
    circle((0.8, -2.4), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "x2")
    content((0.8, -2.4), text(8pt, $x_2$))
    circle((2.2, -2.4), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "C2")
    content((2.2, -2.4), text(8pt, $C_2$))
    
    // Edges
    line("root", "L")
    line("root", "R")
    line("L", "x1")
    line("L", "C1")
    line("R", "x2")
    line("R", "C2")
    
    // Cut line
    line((0, 0.3), (0, -2.7), stroke: (paint: cut-color, dash: "dashed", thickness: 1.5pt))
    content((-1.5, -3), text(8pt, cut-color)[$A = {x_1, C_1}$])
    content((1.5, -3), text(8pt, cut-color)[$overline(A) = {x_2, C_2}$])
  }),
  caption: [Rank decomposition with cut at the root edge.],
)

=== Step 1: Compute Cut Matrix and Rank

The cut matrix $M_A$ for partition $A = {x_1, C_1}$ vs $overline(A) = {x_2, C_2}$:

#v(20pt)
#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    [], [$x_2$], [$C_2$],
    [$x_1$], [0], [1],
    [$C_1$], [1], [0],
  ),
  caption: [Cut matrix $M_A$. Entry is 1 if edge exists in incidence graph.],
)

Over $bb(F)_2$: $op("rank")(M_A) = 2$ (rows are linearly independent). So we have $2^2 = 4$ equivalence classes.

_Note: This small example doesn't show exponential savings. The key is that for $K_(n,n)$, rank stays 1 while cut size grows to $n$._

=== Step 2: DP from Leaves to Root

The DP computes a table $T_v [c]$ at each node $v$, indexed by equivalence class $c$.

*Leaf nodes* (variables $x_1, x_2$): For each assignment $x_i in {0, 1}$, compute its signature.

#v(20pt)
#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([Node], [Assignment], [Table $T$]),
    [$x_1$], [$x_1 = 0$], [$T_(x_1)[sigma_0] = 1$],
    [], [$x_1 = 1$], [$T_(x_1)[sigma_1] = 1$],
    [$x_2$], [$x_2 = 0$], [$T_(x_2)[sigma'_0] = 1$],
    [], [$x_2 = 1$], [$T_(x_2)[sigma'_1] = 1$],
  ),
  caption: [Leaf tables: count 1 for each possible assignment.],
)

*Leaf nodes* (clauses $C_1, C_2$): Clauses contribute constraint information.

*Internal nodes*: Merge child tables by summing over compatible signatures.

=== Step 3: Detailed Merge at Root

At the root, we combine left subtree (covering $x_1, C_1$) and right subtree (covering $x_2, C_2$).

*Left subtree result* $T_L$: counts of $(x_1)$ assignments that satisfy $C_1$'s constraints

#v(20pt)
#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([$x_1$], [Status of $C_1$], [Count]),
    [0], [$C_1$ needs $x_2 = 1$], [1],
    [1], [$C_1$ satisfied], [1],
  ),
  caption: [Left subtree: how $x_1$ affects $C_1 = (x_1 or x_2)$.],
)

*Right subtree result* $T_R$: counts of $(x_2)$ assignments that satisfy $C_2$'s constraints

#v(20pt)
#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([$x_2$], [Status of $C_2$], [Count]),
    [0], [$C_2$ satisfied], [1],
    [1], [$C_2$ needs $x_1 = 0$], [1],
  ),
  caption: [Right subtree: how $x_2$ affects $C_2 = (not x_1 or not x_2)$.],
)

=== Step 4: Final Contraction

Now we contract $T_L$ and $T_R$, checking which combinations satisfy *both* clauses:

#v(20pt)
#figure(
  table(
    columns: 5,
    stroke: 0.5pt,
    table.header([$x_1$], [$x_2$], [$C_1$ ok?], [$C_2$ ok?], [Valid?]),
    [0], [0], [needs $x_2=1$ #sym.times], [#sym.checkmark], [#sym.times],
    [0], [1], [#sym.checkmark], [#sym.checkmark], [#sym.checkmark],
    [1], [0], [#sym.checkmark], [#sym.checkmark], [#sym.checkmark],
    [1], [1], [#sym.checkmark], [needs $x_1=0$ #sym.times], [#sym.times],
  ),
  caption: [Compatibility check: only 2 assignments satisfy both clauses.],
)

*Result*: $\#op("SAT")(phi) = 2$

=== The Tensor Network View

The computation is a tensor contraction:

$ \#op("SAT") = sum_c T_L [c] dot T_R [c] $

#figure(
  canvas(length: 1cm, {
    import draw: *
    
    let tensor-fill = rgb("#a8d5ba")
    
    // Left tensor
    circle((-1.5, 0), radius: 0.5, fill: tensor-fill, stroke: black, name: "L")
    content((-1.5, 0), text(10pt)[$T_L$])
    
    // Bond (equivalence class index)
    line("L", (0.5, 0), stroke: red + 1.5pt, name: "bond")
    content((0, 0.4), text(9pt, red)[bond dim $= 2^r$])
    
    // Right tensor
    circle((1.5, 0), radius: 0.5, fill: tensor-fill, stroke: black, name: "R")
    content((1.5, 0), text(10pt)[$T_R$])
    
    // External legs (summed over)
    line((-1.5, 0.5), (-1.5, 1.2), stroke: gray)
    content((-1.5, 1.5), text(8pt, gray)[$x_1$])
    
    line((1.5, 0.5), (1.5, 1.2), stroke: gray)
    content((1.5, 1.5), text(8pt, gray)[$x_2$])
  }),
  caption: [Tensor network for \#SAT. Bond dimension $= 2^(op("rank")(M_A))$. For $K_(n,n)$, this is $2^1 = 2$ regardless of $n$!],
) <fig:tn-count>

#v(20pt)
*Key insight*: For $K_(n,n)$ with $n$ variables per side, the bond dimension is $2^1 = 2$, not $2^n$. This is the exponential saving from rank-width!

= When Does Rank-Width Help? (Summary)

== Similar Performance

For *sparse* graphs (grid-like, tree-like):
- Treewidth $approx$ Rank-width
- Both approaches give similar complexity

== Exponential Improvement

For *dense* graphs with *low-rank structure*:
- Complete graphs $K_n$: treewidth $= n-1$, rank-width $= 1$
- Complete bipartite $K_(n,n)$: treewidth $= n$, rank-width $= 1$
- Random dense graphs: often have low rank-width

= Connection to Tree Decomposition Framework

#figure(
  table(
    columns: 2,
    stroke: none,
    table.hline(),
    table.header([*Tensor Network Concept*], [*Rank-Width Analogue*]),
    table.hline(),
    [Tree decomposition], [Rank decomposition (binary tree)],
    [Bag (set of vertices)], [Partition induced by tree edge],
    [Separator size], [Rank of cut matrix],
    [Bond dimension], [$2^(op("rank"))$ equivalence classes],
    [Line graph], [Incidence graph (for SAT)],
    table.hline(),
  ),
  caption: [Mapping between tensor network and rank-width concepts.],
) <tab:mapping>

= Practical Implications

For tensor network contraction:

+ *When treewidth is small*: Use standard tree decomposition $arrow.r$ tensor contraction.

+ *When treewidth is large but rank-width is small*: Consider rank-decomposition-based contraction. The effective bond dimension becomes $2^(r w)$ instead of $2^(t w)$.

+ *Computing rank-width*: This is harder than computing treewidth, but approximation algorithms exist @computing-rw.

The main result of @ganian2010 shows that \#SAT can be solved in time $O(2^(3 r w) dot n^2 dot m)$ where $r w$ is the rank-width, which can be exponentially better than clique-width-based methods (since clique-width $<= 2^(r w + 1) - 1$).

= Conclusion

Rank-width provides an alternative to treewidth for parameterizing the complexity of \#SAT and related problems. For tensor network practitioners, the key insight is that rank-width measures the "entanglement" across cuts in terms of $bb(F)_2$-rank rather than raw separator size. This can lead to exponential improvements for graph families where treewidth is large but the adjacency structure has low rank.

#bibliography("refs.bib", title: "References", style: "ieee")

