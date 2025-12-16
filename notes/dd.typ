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

= A Concrete \#SAT Example

Consider the Boolean formula:
$ phi = (x_1 or x_2) and (x_2 or x_3) and (not x_1 or x_3) $

== Step 1: Build the Incidence Graph

The *incidence graph* $I(phi)$ has:
- *Vertices*: Variables ${x_1, x_2, x_3}$ and clauses ${C_1, C_2, C_3}$
- *Edges*: Connect variable to clause if variable appears in clause

#figure(
  canvas(length: 1cm, {
    import draw: *
    
    let var-color = rgb("#4a90d9")
    let clause-color = rgb("#e07b53")
    let d = 2.0
    
    // Variable nodes (circles)
    for (pos, label, nm) in (((0, 0), $x_1$, "x1"), ((d, 0), $x_2$, "x2"), ((d/2, -d), $x_3$, "x3")) {
      circle(pos, radius: 0.35, fill: var-color.lighten(70%), stroke: var-color, name: nm)
      content(pos, text(9pt, label))
    }
    
    // Clause nodes (rectangles)
    for (pos, label, nm) in (((d/2, 0), $C_1$, "C1"), ((d, -d/2), $C_2$, "C2"), ((0, -d/2), $C_3$, "C3")) {
      rect((pos.at(0) - 0.35, pos.at(1) - 0.25), (pos.at(0) + 0.35, pos.at(1) + 0.25), 
           fill: clause-color.lighten(70%), stroke: clause-color, name: nm)
      content(pos, text(9pt, label))
    }
    
    // Edges
    set-style(stroke: gray + 0.8pt)
    line("x1", "C1")
    line("C1", "x2")
    line("x2", "C2")
    line("C2", "x3")
    line("x3", "C3")
    line("C3", "x1")
  }),
  caption: [Incidence graph for $phi = (x_1 or x_2) and (x_2 or x_3) and (not x_1 or x_3)$. Variables (blue circles) connect to clauses (orange rectangles) they appear in.],
) <fig:incidence>

== Step 2: Build the Rank Decomposition

We build a binary tree with leaves corresponding to vertices of $I(phi)$:

#figure(
  canvas(length: 1cm, {
    import draw: *
    
    let node-fill = rgb("#f0f0f0")
    let leaf-fill = rgb("#e8f4e8")
    let cut-color = rgb("#e74c3c")
    let dx = 1.8
    let dy = 1.2
    
    // Level 0 (root)
    circle((0, 0), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "root")
    
    // Level 1
    circle((-dx, -dy), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "L1")
    circle((dx, -dy), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "R1")
    
    // Level 2 - leaves and internal nodes
    circle((-dx - dx/2, -2*dy), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "n-x1")
    content((-dx - dx/2, -2*dy), text(8pt, $x_1$))
    
    circle((-dx + dx/2, -2*dy), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "L2")
    
    circle((dx - dx/2, -2*dy), radius: 0.2, fill: node-fill, stroke: black + 0.6pt, name: "R2")
    
    circle((dx + dx/2, -2*dy), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "n-x3")
    content((dx + dx/2, -2*dy), text(8pt, $x_3$))
    
    // Level 3 - all leaves
    circle((-dx + dx/2 - 0.6, -3*dy), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "n-C1")
    content((-dx + dx/2 - 0.6, -3*dy), text(8pt, $C_1$))
    
    circle((-dx + dx/2 + 0.6, -3*dy), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "n-C3")
    content((-dx + dx/2 + 0.6, -3*dy), text(8pt, $C_3$))
    
    circle((dx - dx/2 - 0.6, -3*dy), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "n-x2")
    content((dx - dx/2 - 0.6, -3*dy), text(8pt, $x_2$))
    
    circle((dx - dx/2 + 0.6, -3*dy), radius: 0.3, fill: leaf-fill, stroke: black + 0.6pt, name: "n-C2")
    content((dx - dx/2 + 0.6, -3*dy), text(8pt, $C_2$))
    
    // Edges
    set-style(stroke: black + 0.6pt)
    line("root", "L1")
    line("root", "R1")
    line("L1", "n-x1")
    line("L1", "L2")
    line("R1", "R2")
    line("R1", "n-x3")
    line("L2", "n-C1")
    line("L2", "n-C3")
    line("R2", "n-x2")
    line("R2", "n-C2")
    
    // Cut line (dashed red)
    line((0, 0.5), (0, -3.8*dy), stroke: (paint: cut-color, dash: "dashed", thickness: 1.5pt))
    
    // Labels for partitions
    content((-dx, -4.2*dy), text(8pt, cut-color)[$A = {x_1, C_1, C_3}$])
    content((dx, -4.2*dy), text(8pt, cut-color)[$overline(A) = {x_2, C_2, x_3}$])
  }),
  caption: [A rank decomposition for the incidence graph. Leaves (green) correspond to vertices. The dashed red line shows the main cut partitioning into $A$ and $overline(A)$.],
) <fig:rank-decomp>

== Step 3: Compute Cut Matrices

For the main cut $A = {x_1, C_1, C_3}$ vs $overline(A) = {x_2, C_2, x_3}$:

#v(20pt)
#figure(
  table(
    columns: 4,
    stroke: 0.5pt,
    [], [$x_2$], [$C_2$], [$x_3$],
    [$x_1$], [0], [0], [0],
    [$C_1$], [1], [0], [0],
    [$C_3$], [0], [0], [1],
  ),
  caption: [Cut matrix $M_A$ for partition $A = {x_1, C_1, C_3}$.],
) <tab:cut-matrix>

The $bb(F)_2$-rank of this matrix is *2* (two linearly independent rows over the binary field).

== Step 4: The DP State Space

Here lies the key insight connecting to tensor networks:

*Treewidth-based algorithms*: The DP state tracks *all possible assignments* to variables in the separator. If the separator has $k$ variables, we have $2^k$ states.

*Rank-width-based algorithms*: The DP state tracks *equivalence classes* of partial assignments. Two partial assignments are equivalent if they produce the same "signature" on the boundary.

=== Understanding Equivalence Classes

Consider the cut $A = {x_1, C_1, C_3}$ vs $overline(A) = {x_2, C_2, x_3}$. The edges crossing this cut are:
- $(C_1, x_2)$: clause $C_1 = (x_1 or x_2)$ connects to $x_2$
- $(C_3, x_3)$: clause $C_3 = (not x_1 or x_3)$ connects to $x_3$

For a partial assignment to $x_1$ (the only variable in $A$), we ask: _what information does the other side need?_

#v(20pt)
#figure(
  table(
    columns: 4,
    stroke: 0.5pt,
    table.header([$x_1$], [Status of $C_1$], [Status of $C_3$], [*Signature*]),
    [0], [needs $x_2 = 1$], [satisfied], [$(0, 1)$],
    [1], [satisfied], [needs $x_3 = 1$], [$(1, 0)$],
  ),
  caption: [The "signature" encodes which boundary constraints remain.],
) <tab:equiv>

The *signature* is a vector in $bb(F)_2^r$ where $r = op("rank")(M_A) = 2$. It encodes:
- For each crossing edge, whether the clause on side $A$ is already satisfied (1) or still needs help from $overline(A)$ (0).

Two partial assignments are *equivalent* if they produce the same signature. In this example, $x_1 = 0$ and $x_1 = 1$ produce different signatures, so we have 2 equivalence classes.

=== SVD Perspective (Tensor Network View)

For tensor network practitioners, the rank decomposition is best understood via SVD. Consider the cut matrix:
$ M_A = mat(0, 0, 0; 1, 0, 0; 0, 0, 1) $
with rows ${x_1, C_1, C_3}$ and columns ${x_2, C_2, x_3}$.

*Over $bb(R)$*, we can compute the SVD: $M_A = U Sigma V^top$ where $op("rank")(M_A) = 2$. This factorizes as:
$ M_A = underbrace(mat(0, 0; 1, 0; 0, 1), U') dot underbrace(mat(1, 0, 0; 0, 0, 1), V'^top) $

*The key insight*: Any vector $bold(v) in bb(F)_2^(|A|)$ indexing rows of $M_A$ produces a "boundary effect" $bold(v)^top M_A in bb(F)_2^(|overline(A)|)$. But since $op("rank")(M_A) = r$, this effect factors through an $r$-dimensional space:
$ bold(v)^top M_A = (bold(v)^top U') dot V'^top $

The *equivalence class* of $bold(v)$ is determined by $bold(v)^top U' in bb(F)_2^r$ — a vector of dimension $r$, not $|A|}$.

#figure(
  canvas(length: 1cm, {
    import draw: *
    
    let box-fill = rgb("#f8f8f8")
    
    // Left: high-dimensional space
    rect((-0.5, -1), (1.5, 1), fill: box-fill, stroke: gray)
    content((0.5, 1.3), text(9pt)[$bb(F)_2^(|A|)$])
    content((0.5, 0), text(8pt)[all partial\ assignments])
    
    // Arrow to middle
    line((1.7, 0), (2.8, 0), mark: (end: "stealth"), stroke: black + 0.8pt)
    content((2.25, 0.4), text(8pt)[$times U'$])
    
    // Middle: low-dimensional space
    rect((3, -0.7), (4.5, 0.7), fill: rgb("#e8f4e8"), stroke: rgb("#4a4"))
    content((3.75, 1.0), text(9pt)[$bb(F)_2^r$])
    content((3.75, 0), text(8pt)[equiv.\ classes])
    
    // Arrow to right
    line((4.7, 0), (5.8, 0), mark: (end: "stealth"), stroke: black + 0.8pt)
    content((5.25, 0.4), text(8pt)[$times V'^top$])
    
    // Right: boundary effect
    rect((6, -1), (8, 1), fill: box-fill, stroke: gray)
    content((7, 1.3), text(9pt)[$bb(F)_2^(|overline(A)|)$])
    content((7, 0), text(8pt)[boundary\ effects])
  }),
  caption: [SVD factorization: partial assignments map to equivalence classes (dimension $r$), which then determine boundary effects. The bottleneck has dimension $r = op("rank")(M_A)$.],
) <fig:svd>

#v(20pt)

The number of distinct equivalence classes is at most $2^(op("rank")(M_A))$, not $2^(|A|)$.

== Step 5: The Tensor Network View

Think of each clause as a tensor:
- $C_1 = (x_1 or x_2)$: A tensor $T_(C_1) [x_1, x_2] in {0, 1}$
- $C_2 = (x_2 or x_3)$: A tensor $T_(C_2) [x_2, x_3] in {0, 1}$
- $C_3 = (not x_1 or x_3)$: A tensor $T_(C_3) [x_1, x_3] in {0, 1}$

The \#SAT problem becomes a *tensor network contraction*:
$ \#op("SAT")(phi) = sum_(x_1, x_2, x_3 in {0,1}) T_(C_1) [x_1, x_2] dot T_(C_2) [x_2, x_3] dot T_(C_3) [x_1, x_3] $

== The Rank-Width Advantage

When contracting across a cut, we don't need to track all $2^k$ boundary configurations. We only need $2^(r w)$ *equivalence classes*.

Over $bb(F)_2$, the equivalence is determined by *linear combinations* of edge indicators. If the cut matrix has rank $r$, there are only $2^r$ distinct "boundary behaviors."

= A Better Example: Dense Formulas

The previous example had treewidth $approx$ rank-width $approx 2$, showing no advantage. Let's see a case where rank-width is *exponentially smaller*.

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

== Why Does This Work?

For any cut in $K_(n,n)$, consider the cut matrix. Since all entries are 1, any partial assignment produces a boundary effect that is either:
- The all-zeros vector (if an even number of edges cross), or
- The all-ones vector (if an odd number of edges cross)

There are only *2 equivalence classes* over $bb(F)_2$, regardless of how many vertices are on each side.

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

