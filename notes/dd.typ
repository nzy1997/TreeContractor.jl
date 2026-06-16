#import "@preview/clear-iclr:0.7.0": *
#import "@preview/cetz:0.4.2": canvas, draw
#import "@preview/ctheorems:1.1.3": thmbox, thmrules

// Theorem environments
#let definition = thmbox("definition", "Definition", fill: rgb("#f0f8ff"), stroke: rgb("#4a90d9"))
#let proposition = thmbox("proposition", "Proposition", fill: rgb("#fff8e8"), stroke: rgb("#c90"))
#let remark = thmbox("remark", "Remark", fill: rgb("#e8f8e8"), stroke: rgb("#4a4"))
#let algorithm = thmbox("algorithm", "Algorithm", fill: rgb("#f5f5f5"), stroke: rgb("#999"))

#show: thmrules

#show: iclr.with(
  title: [A Tutorial on Rank-Width for \#SAT],
  authors: (
    (
      names: ([Anonymous],),
      affilation: [Tutorial Notes],
      address: [],
      email: "cacate0129@gmail.com",
    ),
  ),
  keywords: ("rank-width", "SAT", "treewidth", "parameterized algorithms"),
  abstract: [
    We introduce rank-width as a graph parameter for the model counting problem \#SAT. The core idea is neighborhood equivalence: partial assignments with the same $bb(F)_2$ characteristic interact identically across any cut, enabling state-space compression from $2^n$ to $2^(op("rw"))$.
  ],
  accepted: none,
  bibliography: none,
)

= Introduction

The *model counting problem* \#SAT asks: given a Boolean formula $phi$, how many satisfying truth assignments exist? While \#P-complete in general, efficient algorithms exist when the formula's structure is restricted.
Ganian, Hliněný, and Obdržálek @ganian2010 provide a \#SAT algorithm with runtime single-exponential in the _rank-width_ of the incidence graph. *Treewidth* measures how "tree-like" a graph is; it can be large for dense graphs, but rank-width can be exponentially smaller for certain graph families.

= Rank Decomposition and Rank-Width

A *rank decomposition* $(T, L)$ of graph $G = (V, E)$ consists of a subcubic tree $T$ (every node has degree $<= 3$) and a bijection $L$ from leaves to $V$. Each edge $e in T$ induces a partition $(A_e, B_e)$ of $V$: removing $e$ splits the tree into two subtrees, with $A_e$ and $B_e$ being the leaves in each.

#figure(
  canvas(length: 1cm, {
    import draw: *
    let leaf-col = rgb("#4a90d9")
    let internal-col = rgb("#888")
    let cut-col = rgb("#e07b53")
    // Tree structure
    circle((0, 0), radius: 0.2, fill: internal-col.lighten(50%), stroke: internal-col, name: "root")
    circle((-2, -1.5), radius: 0.2, fill: internal-col.lighten(50%), stroke: internal-col, name: "i1")
    circle((2, -1.5), radius: 0.2, fill: internal-col.lighten(50%), stroke: internal-col, name: "i2")
    // Leaves
    circle((-3, -3), radius: 0.3, fill: leaf-col.lighten(70%), stroke: leaf-col, name: "l1")
    content((-3, -3), text(7pt)[$v_1$])
    circle((-1, -3), radius: 0.3, fill: leaf-col.lighten(70%), stroke: leaf-col, name: "l2")
    content((-1, -3), text(7pt)[$v_2$])
    circle((1, -3), radius: 0.3, fill: leaf-col.lighten(70%), stroke: leaf-col, name: "l3")
    content((1, -3), text(7pt)[$v_3$])
    circle((3, -3), radius: 0.3, fill: leaf-col.lighten(70%), stroke: leaf-col, name: "l4")
    content((3, -3), text(7pt)[$v_4$])
    // Edges
    line("root", "i1")
    line("root", "i2")
    line("i1", "l1")
    line("i1", "l2")
    line("i2", "l3")
    line("i2", "l4")
    // Cut edge highlight
    line((-0.15, -0.15), (0.15, 0.15), stroke: cut-col + 2pt)
    content((0.8, 0.3), text(8pt, cut-col)[cut edge $e$])
    // Partition labels
    content((-2, -3.8), text(8pt)[$A_e = {v_1, v_2}$])
    content((2, -3.8), text(8pt)[$B_e = {v_3, v_4}$])
  }),
  // caption: [A rank decomposition. Removing edge $e$ partitions the leaves into $A_e$ and $B_e$.],
) <fig:rank-decomp>

#v(10pt)
For partition $(A, B)$, the *cut matrix* $M_A in bb(F)_2^(|A| times |B|)$ has $M_A [u,v] = 1$ iff $(u,v) in E$. The *rank-width* is:
$ op("rw")(G) = min_((T,L)) max_(e in T) op("rank")_(bb(F)_2)(M_(A_e)) $

= Neighborhood Equivalence

Let $G = (V, E)$ be a graph. Given a partition $(A, B)$ of $V$ and a subset $Y subset.eq A$, the *neighborhood* of $Y$ relative to $B$ is defined as:
$ N_B (Y) = {v in B : exists y in Y, {y, v} in E}. $

Two subsets $Y, Y' subset.eq A$ are said to be *neighborhood-equivalent* with respect to the partition $(A, B)$ if and only if $N_B (Y) = N_B (Y')$.

#figure(
  canvas(length: 1cm, {
    import draw: *
    let inside-col = rgb("#4a90d9")
    let boundary-col = rgb("#e07b53")
    let class1-col = rgb("#9b59b6")  // purple for class 1
    let class2-col = rgb("#27ae60")  // green for class 2
    // Inside region
    rect((-4, -2.5), (1, 2.5), fill: inside-col.lighten(90%), stroke: inside-col)
    content((-1.5, 2.8), text(9pt, inside-col)[$A$])
    // Boundary region
    rect((1.5, -2.5), (4.5, 2.5), fill: boundary-col.lighten(90%), stroke: boundary-col)
    content((3, 2.8), text(9pt, boundary-col)[$B$])
    // Vertices inside - 4 vertices
    circle((-3, 1.5), radius: 0.3, fill: class1-col.lighten(50%), stroke: class1-col + 1.5pt, name: "a1")
    content((-3, 1.5), text(8pt)[$a_1$])
    circle((0, 2), radius: 0.3, fill: class1-col.lighten(50%), stroke: class1-col + 1.5pt, name: "a2")
    content((0, 2), text(8pt)[$a_2$])
    circle((-3, 0), radius: 0.3, fill: class1-col.lighten(50%), stroke: class1-col + 1.5pt, name: "a3")
    content((-3, 0), text(8pt)[$a_3$])
    circle((-0.5, -0.5), radius: 0.3, fill: class2-col.lighten(50%), stroke: class2-col + 1.5pt, name: "a4")
    content((-0.5, -0.5), text(8pt)[$a_4$])
    // Vertices boundary - 2 vertices
    circle((3, 1), radius: 0.3, fill: boundary-col.lighten(50%), stroke: boundary-col, name: "b1")
    content((3, 1), text(8pt)[$b_1$])
    circle((3, -0.5), radius: 0.3, fill: boundary-col.lighten(50%), stroke: boundary-col, name: "b2")
    content((3, -0.5), text(8pt)[$b_2$])
    // Cross edges: a1,a2,a3 -> b1 only; a4 -> b1,b2
    line("a1", "b1", stroke: class1-col + 1pt)
    line("a2", "b1", stroke: class1-col + 1pt)
    line("a3", "b1", stroke: class1-col + 1pt)
    line("a4", "b1", stroke: class2-col + 1pt)
    line("a4", "b2", stroke: class2-col + 1pt)
    // Annotation
    content((-1.5, -1.8), text(8pt, class1-col)[Class 1: $N_B = {b_1}$])
    content((3, -1.8), text(8pt, class2-col)[Class 2: $N_B = {b_1, b_2}$])
  }),
) <fig:cut-boundary>
#v(10pt)

In @fig:cut-boundary, vertices $a_1, a_2, a_3$ (purple) all connect only to $b_1$, so they share the same neighborhood $N_B = {b_1}$ and belong to *Class 1*. Vertex $a_4$ (green) connects to both $b_1$ and $b_2$, so $N_B ({a_4}) = {b_1, b_2}$ --- *Class 2*.

#v(10pt)
#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([*Neighborhood*], [*Subsets*], [*Count*]),
    [$emptyset$], [$emptyset$], [1],
    [${b_1}$], [${a_1}, {a_2}, {a_3}, {a_1, a_2}, {a_1, a_3}, {a_2, a_3}, {a_1, a_2, a_3}$], [7],
    [${b_1, b_2}$], [any subset containing $a_4$], [8],
  ),
)

This property enables efficient \#SAT computation. In a dynamic programming approach, at each step we track how partial assignments interact with clauses not yet fully determined. Partial assignments with identical neighborhoods combine identically with any assignment to the remaining variables. Consequently, instead of tracking $2^(|A_e|)$ partial assignments separately, we track only $O(2^k)$ equivalence classes, where $k$ is the rank-width.

== Example: $K_(3,3)$ 

To apply rank-width to \#SAT, we construct the *incidence graph* of the CNF formula: vertices represent variables and clauses, with an edge between variable $x_i$ and clause $C_j$ if $x_i$ appears (positively or negatively) in $C_j$. The rank-width of this incidence graph determines the algorithm's efficiency.

Consider a CNF formula where *every variable appears in every clause*:
$ phi = C_1 and C_2 and C_3 $
where:
- $C_1 = (x_1 or x_2 or x_3)$
- $C_2 = (not x_1 or x_2 or x_3)$
- $C_3 = (x_1 or not x_2 or not x_3)$

With variables $X = {x_1, x_2, x_3}$ and clauses $C = {C_1, C_2, C_3}$, each variable appears (positively or negatively) in every clause, so the incidence graph is $K_(3,3)$.

#figure(
  canvas(length: 1cm, {
    import draw: *
    let var-col = rgb("#4a90d9")
    let clause-col = rgb("#e07b53")
    // Variables (left side)
    circle((-2, 1.5), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x1")
    content((-2, 1.5), text(9pt)[$x_1$])
    circle((-2, 0), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x2")
    content((-2, 0), text(9pt)[$x_2$])
    circle((-2, -1.5), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x3")
    content((-2, -1.5), text(9pt)[$x_3$])
    // Clauses (right side)
    circle((2, 1.5), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c1")
    content((2, 1.5), text(9pt)[$C_1$])
    circle((2, 0), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c2")
    content((2, 0), text(9pt)[$C_2$])
    circle((2, -1.5), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c3")
    content((2, -1.5), text(9pt)[$C_3$])
    // All edges (complete bipartite)
    for v in ("x1", "x2", "x3") {
      for c in ("c1", "c2", "c3") {
        line(v, c, stroke: gray + 0.8pt)
      }
    }
    // Labels
    content((-2, 2.3), text(8pt, var-col)[Variables])
    content((2, 2.3), text(8pt, clause-col)[Clauses])
  }),
  caption: [Incidence graph $K_(3,3)$: every variable connected to every clause.],
) <fig:k33-graph>

#v(10pt)
The cut matrix separating variables $X$ from clauses $C$ is:
$ M = mat(1, 1, 1; 1, 1, 1; 1, 1, 1) $

Over $bb(F)_2$, this matrix has *rank 1* (all rows identical). The equivalence class of a truth assignment is determined by its *neighborhood characteristic*.

=== The Neighborhood Characteristic

Recall that our goal is to count satisfying assignments by dynamic programming on the rank decomposition. At each node, the algorithm processes a subset $A subset.eq V$ of vertices (variables and clauses) and must track how partial assignments to $A$ interact with the unprocessed boundary $B = V - A$.

A naive approach would enumerate all $2^(|A|)$ possible subsets of TRUE variables, leading to exponential state space. However, observe that the subsequent computation depends only on *how* the partial assignment connects to $B$, not on *which* specific variables are selected. Two partial assignments that induce identical adjacency patterns with $B$ will combine identically with any completion on $B$. This observation motivates the following definition.

#definition("Neighborhood Characteristic")[
  Given a subset $S subset.eq A$ with characteristic vector $bold(v) in bb(F)_2^n$ (where $v_i = 1$ iff $x_i in S$), the *neighborhood characteristic* of $S$ is:
  $ bold(n)_S = bold(v)^top M in bb(F)_2^m $
  where $M$ is the cut matrix of the partition $(A, B)$.
]

The $j$-th entry $(bold(n)_S)_j = plus.big_(i : x_i in S) M[i,j]$ counts (mod 2) the number of vertices in $S$ adjacent to the $j$-th boundary vertex.

The choice of the binary field $bb(F)_2$ is essential for compositionality. When joining partial assignments $S_1$ and $S_2$ from disjoint subtrees, their combined characteristic satisfies:
$ bold(n)_(S_1 union S_2) = bold(n)_(S_1) xor bold(n)_(S_2) $
This closure property ensures that $bb(F)_2$-equivalent partial assignments remain equivalent after composition: if $bold(n)_(S_1) = bold(n)_(S'_1)$, then $bold(n)_(S_1 union S_2) = bold(n)_(S'_1 union S_2)$ for any $S_2$.

Two subsets are $bb(F)_2$-equivalent if and only if they have identical neighborhood characteristics:
$ S tilde.eq_(bb(F)_2) S' <==> bold(n)_S = bold(n)_(S') <==> bold(v) - bold(v)' in ker(M) $
The number of equivalence classes equals $2^(op("rank")(M))$ @bui-xuan2011. Consequently, the DP algorithm tracks at most $2^(op("rw")(G))$ equivalence classes per node, rather than $2^(|A|)$ individual subsets.

*Application to \#SAT.* When $S$ represents the set of TRUE variables and $B$ the set of clauses, $(bold(n)_S)_j$ counts (mod 2) how many TRUE variables appear in clause $C_j$. This encodes *incidence* (which variables appear in which clauses) rather than *satisfaction* (whether clauses evaluate to true). Literal polarity---whether a variable appears as $x_i$ or $not x_i$---is tracked separately in the DP state; the neighborhood characteristic captures interaction patterns across the cut, while satisfaction verification requires additional bookkeeping described in subsequent sections.

=== Computing the Neighborhood Characteristic

Let $bold(v) in bb(F)_2^3$ be the characteristic vector of TRUE variables (e.g., $bold(v) = (1,1,0)^top$ means $x_1 = x_2 = "TRUE", x_3 = "FALSE"$).

For our dense formula with $M = J$ (all-ones matrix):
$ bold(n)_S = bold(v)^top M = bold(v)^top mat(1, 1, 1; 1, 1, 1; 1, 1, 1) = (v_1 xor v_2 xor v_3) dot (1, 1, 1)^top $

Since every variable appears in every clause, the neighborhood characteristic depends only on the *parity* of $|S|$: odd-sized sets yield $(1,1,1)^top$, even-sized sets yield $(0,0,0)^top$.

#v(10pt)
#figure(
  table(
    columns: 5,
    stroke: 0.5pt,
    table.header([*TRUE vars $S$*], [*$bold(v)^top$*], [*$|S| mod 2$*], [*Neighborhood $bold(n)_S^top$*], [*Class*]),
    [$emptyset$], [$(0,0,0)$], [$0$], [$(0,0,0)$], [0],
    [${x_1}$], [$(1,0,0)$], [$1$], [$(1,1,1)$], [1],
    [${x_2}$], [$(0,1,0)$], [$1$], [$(1,1,1)$], [1],
    [${x_3}$], [$(0,0,1)$], [$1$], [$(1,1,1)$], [1],
    [${x_1, x_2}$], [$(1,1,0)$], [$0$], [$(0,0,0)$], [0],
    [${x_1, x_3}$], [$(1,0,1)$], [$0$], [$(0,0,0)$], [0],
    [${x_2, x_3}$], [$(0,1,1)$], [$0$], [$(0,0,0)$], [0],
    [${x_1, x_2, x_3}$], [$(1,1,1)$], [$1$], [$(1,1,1)$], [1],
  ),
  caption: [The neighborhood characteristic $bold(n)_S$ determines the equivalence class.],
) <tab:k33-equiv>

#remark[
  The $bb(F)_2$ equivalence is coarser than set-theoretic neighborhood equivalence. For instance, $emptyset$ has $N_B = emptyset$ while ${x_1, x_2}$ has $N_B = {C_1, C_2, C_3}$; these are not set-theoretically equivalent, yet both have characteristic $(0,0,0)^top$ since ${x_1, x_2}$ contributes $1 xor 1 = 0$ to each clause. This coarser grouping suffices because the DP algorithm only requires that partial assignments with identical characteristics combine identically with any boundary assignment.
]

== Counterexample: Sparse Formula --- No Compression

Now consider a CNF formula where *each variable appears in exactly one clause*:
$ psi = C_1 and C_2 and C_3 $
where $C_1 = (x_1), C_2 = (x_2), C_3 = (x_3)$.

With variables $X = {x_1, x_2, x_3}$ and clauses $C = {C_1, C_2, C_3}$, each variable appears in exactly one clause, so the incidence graph is a *perfect matching* $M_3$.

#figure(
  canvas(length: 1cm, {
    import draw: *
    let var-col = rgb("#4a90d9")
    let clause-col = rgb("#e07b53")
    // Variables (left side)
    circle((-2, 1.5), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x1")
    content((-2, 1.5), text(9pt)[$x_1$])
    circle((-2, 0), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x2")
    content((-2, 0), text(9pt)[$x_2$])
    circle((-2, -1.5), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x3")
    content((-2, -1.5), text(9pt)[$x_3$])
    // Clauses (right side)
    circle((2, 1.5), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c1")
    content((2, 1.5), text(9pt)[$C_1$])
    circle((2, 0), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c2")
    content((2, 0), text(9pt)[$C_2$])
    circle((2, -1.5), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c3")
    content((2, -1.5), text(9pt)[$C_3$])
    // Only matching edges (x_i to C_i)
    line("x1", "c1", stroke: gray + 1.2pt)
    line("x2", "c2", stroke: gray + 1.2pt)
    line("x3", "c3", stroke: gray + 1.2pt)
    // Labels
    content((-2, 2.3), text(8pt, var-col)[Variables])
    content((2, 2.3), text(8pt, clause-col)[Clauses])
  }),
  caption: [Incidence graph $M_3$: each variable connected to exactly one clause.],
) <fig:m3-graph>
#v(10pt)
The cut matrix is the *identity*: $M = I_3$ with $op("rank") = 3$ (full rank). Since $bold(n)_S = bold(v)^top I = bold(v)$, every assignment has a unique neighborhood.
#v(10pt)
#figure(
  table(
    columns: 4,
    stroke: 0.5pt,
    table.header([*TRUE variables $S$*], [*$bold(v)^top$*], [*Neighborhood $bold(n)_S^top$*], [*Class*]),
    [$emptyset$], [$(0,0,0)$], [$(0,0,0)$], [0],
    [${x_1}$], [$(1,0,0)$], [$(1,0,0)$], [1],
    [${x_2}$], [$(0,1,0)$], [$(0,1,0)$], [2],
    [${x_3}$], [$(0,0,1)$], [$(0,0,1)$], [3],
    [${x_1, x_2}$], [$(1,1,0)$], [$(1,1,0)$], [4],
    [${x_1, x_3}$], [$(1,0,1)$], [$(1,0,1)$], [5],
    [${x_2, x_3}$], [$(0,1,1)$], [$(0,1,1)$], [6],
    [${x_1, x_2, x_3}$], [$(1,1,1)$], [$(1,1,1)$], [7],
  ),
  caption: [Sparse formula: each assignment has a unique neighborhood --- no compression.],
) <tab:matching-equiv>

= The DP Algorithm
#text(fill: blue)[HM: below not polished yet.]

The Ganian et al. algorithm @ganian2010 processes the rank decomposition tree bottom-up, using neighborhood equivalence to compress the state space.

== State Representation

The actual algorithm @ganian2010 uses *signed graphs* to distinguish positive and negative literal occurrences. The signed graph $F_phi$ has two edge sets: $E^+$ (variable appears positively in clause) and $E^-$ (variable appears negatively). Crucially, *both signed components share the same decomposition tree* --- they differ only in their labeling functions at each node.

Before describing the full state structure, we must introduce *orthogonality*, which encodes partial clause satisfaction. Given vectors $bold(u), bold(v) in bb(F)_2^t$, they are *orthogonal* (written $bold(u) perp bold(v)$) if $bold(u)^top bold(v) = 0$ in $bb(F)_2$. In the DP algorithm, a clause's label vector represents which literal types (positive/negative) can satisfy it. A subspace $Pi^+$ (expected positive contributions) is sufficient to satisfy that clause iff the clause label is *not orthogonal* to $Pi^+$. If a clause label is orthogonal to both $Pi^+$ and $Pi^-$ (expected positive and negative contributions), then no satisfying literal can exist for that clause in this partial assignment, making the state invalid.

At each node $z$ of this shared parse tree:

- *Labeling*: Each variable or clause vertex $v$ has two label vectors $lambda^+(v), lambda^-(v) in bb(F)_2^t$. Variable labels come from the columns of the signed cut matrices (M^+ and M^-). Clause labels are chosen as a basis for the row space of these matrices, encoding which variable combinations can satisfy each clause.

- *State*: A 4-tuple of subspaces $(Sigma^+, Sigma^-, Pi^+, Pi^-)$ over $bb(F)_2^t$ where $t = max(t^+, t^-)$ is the signed rank-width.
  - $Sigma^+$: subspace generated by $lambda^+(v)$ for all variables $v$ assigned TRUE
  - $Sigma^-$: subspace generated by $lambda^-(v)$ for all variables $v$ assigned FALSE
  - $Pi^+, Pi^-$: *unsatisfied clause coverage requirements* --- which clause label patterns still need to be covered by the remaining variables to satisfy clauses
- *Table*: $T_z [Sigma^+, Sigma^-, Pi^+, Pi^-]$ = count of partial assignments of *shape* $(Sigma^+, Sigma^-, Pi^+, Pi^-)$

*Shape Definition*: An assignment $nu$ has shape $(Sigma^+, Sigma^-, Pi^+, Pi^-)$ iff:
+ $Sigma^+$ and $Sigma^-$ are the exact subspaces spanned by the labels of TRUE and FALSE variables respectively
+ Every clause $c$ in the processed subtree satisfies: either $lambda(c) not perp Sigma^+$ (satisfied by positive literals), or $lambda(c) not perp Sigma^-$ (satisfied by negative literals), or $lambda(c) not perp Pi^+$ or $lambda(c) not perp Pi^-$ (expected to be satisfied by remaining variables)

=== The Expectation Mechanism

Tracking *expectations* allows efficient composition. When joining two subtrees, it is unnecessary to enumerate all ways clauses might be satisfied; one need only verify that the expectations from one side match the contributions from the other.

#proposition("Expectation Matching")[
  An assignment $nu$ is satisfying for $phi$ iff there exist subspaces such that the left partial assignment has shape $(Sigma^+, Sigma^-, Pi^+, Pi^-)$ and the right has shape $(Pi^+, Pi^-, Sigma^+, Sigma^-)$ --- the expectations and contributions are swapped.
]

This "expectation" technique, introduced in @bui-xuan2011 for dominating set, is essential for achieving single-exponential runtime in rank-width.

== Processing Nodes

The algorithm processes three types of nodes:

=== Leaf Nodes

The algorithm processes *two types* of leaves differently:

*Variable leaf* $ell$: For each variable, we record both possible assignments ($ell = 0$ or $ell = 1$). Each assignment contributes to a neighborhood class based on its label vector. Writing $op("span")(bold(1))$ for the 1-dimensional subspace spanned by the label vector:
- $T[op("span")(bold(1)), emptyset, Pi^+, Pi^-]$ += 1 for all $Pi^+, Pi^-$ (assignment $ell = 1$)
- $T[emptyset, op("span")(bold(1)), Pi^+, Pi^-]$ += 1 for all $Pi^+, Pi^-$ (assignment $ell = 0$)

*Clause leaf* $c$: A clause starts unsatisfied and requires either a satisfying literal or an "expectation" that it will be satisfied later:
- $T[emptyset, emptyset, Pi^+, Pi^-]$ += 1 only for subspaces $Pi^+, Pi^-$ where at least one is *not orthogonal* to the label vector $bold(1)$

The orthogonality condition encodes clause satisfaction: if a clause's label is orthogonal to all expected satisfying contributions, the clause cannot be satisfied.

=== Internal Nodes (Join)

For an internal node with children corresponding to edges $e_1, e_2$:
$ T_e [bold(s)] = sum_(bold(s)_1 xor bold(s)_2 = bold(s)) T_(e_1)[bold(s)_1] dot T_(e_2)[bold(s)_2] $

This is a *convolution over $bb(F)_2$* --- the XOR reflects that neighborhoods combine additively in $bb(F)_2$.

=== Root Node (Final Count)

At the root, sum over all states where *no expectations remain* (all clauses must be satisfied without relying on future contributions):
$ \#op("SAT")(phi) = sum_(Sigma^+, Sigma^-) T_("root")[Sigma^+, Sigma^-, emptyset, emptyset] $

The condition $Pi^+ = Pi^- = emptyset$ ensures that every clause has been satisfied by some literal in the complete assignment --- no clause is left "expecting" satisfaction from a non-existent remaining part of the formula.

// == Worked Example: Full DP on $K_(3,3)$ with Signed Graphs

// We demonstrate the complete algorithm on the $K_(3,3)$ formula:
// $ phi = (x_1 or x_2 or x_3) and (not x_1 or x_2 or x_3) and (x_1 or not x_2 or not x_3) $

// This formula has exactly 5 satisfying assignments: ${x_2}, {x_3}, {x_1, x_2}, {x_1, x_3}, {x_1, x_2, x_3}$.

// === Step 1: Construct the Signed Graph

// The signed incidence graph $F_phi$ has vertices $V = {x_1, x_2, x_3, C_1, C_2, C_3}$ and two edge sets based on literal polarity:

// - $E^+ = {(x_1, C_1), (x_2, C_1), (x_3, C_1), (x_2, C_2), (x_3, C_2), (x_1, C_3)}$
// - $E^- = {(x_1, C_2), (x_2, C_3), (x_3, C_3)}$

// #figure(
//   canvas(length: 1cm, {
//     import draw: *
//     let var-col = rgb("#4a90d9")
//     let clause-col = rgb("#e07b53")
//     let pos-col = rgb("#27ae60")
//     let neg-col = rgb("#9b59b6")
//     // Variables (left side)
//     circle((-2.5, 1.5), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x1")
//     content((-2.5, 1.5), text(9pt)[$x_1$])
//     circle((-2.5, 0), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x2")
//     content((-2.5, 0), text(9pt)[$x_2$])
//     circle((-2.5, -1.5), radius: 0.35, fill: var-col.lighten(70%), stroke: var-col + 1.5pt, name: "x3")
//     content((-2.5, -1.5), text(9pt)[$x_3$])
//     // Clauses (right side)
//     circle((2.5, 1.5), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c1")
//     content((2.5, 1.5), text(9pt)[$C_1$])
//     circle((2.5, 0), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c2")
//     content((2.5, 0), text(9pt)[$C_2$])
//     circle((2.5, -1.5), radius: 0.35, fill: clause-col.lighten(70%), stroke: clause-col + 1.5pt, name: "c3")
//     content((2.5, -1.5), text(9pt)[$C_3$])
//     // Positive edges (solid green) - C1: all positive
//     line("x1", "c1", stroke: pos-col + 1pt)
//     line("x2", "c1", stroke: pos-col + 1pt)
//     line("x3", "c1", stroke: pos-col + 1pt)
//     // C2: x2, x3 positive
//     line("x2", "c2", stroke: pos-col + 1pt)
//     line("x3", "c2", stroke: pos-col + 1pt)
//     // C3: x1 positive
//     line("x1", "c3", stroke: pos-col + 1pt)
//     // Negative edges (dashed purple) - C2: x1 negative
//     line("x1", "c2", stroke: (paint: neg-col, thickness: 1pt, dash: "dashed"))
//     // C3: x2, x3 negative
//     line("x2", "c3", stroke: (paint: neg-col, thickness: 1pt, dash: "dashed"))
//     line("x3", "c3", stroke: (paint: neg-col, thickness: 1pt, dash: "dashed"))
//     // Legend
//     content((0, -2.5), text(8pt)[#text(pos-col)[solid] = $E^+$, #text(neg-col)[dashed] = $E^-$])
//   }),
//   caption: [Signed incidence graph for $K_(3,3)$: positive and negative literal edges.],
// ) <fig:k33-signed>

// === Step 2: Compute Signed Cut Matrices

// For the partition $(X, C) = ({x_1, x_2, x_3}, {C_1, C_2, C_3})$, the signed cut matrices are:

// $ M^+ = mat(1, 0, 1; 1, 1, 0; 1, 1, 0) quad M^- = mat(0, 1, 0; 0, 0, 1; 0, 0, 1) $

// where $M^+[i,j] = 1$ iff $x_i$ appears positively in $C_j$, and $M^-[i,j] = 1$ iff $x_i$ appears negatively in $C_j$.

// Over $bb(F)_2$: rows 2 and 3 of $M^+$ are identical, so $op("rank")(M^+) = 2$. Similarly, $op("rank")(M^-) = 2$. Thus $t = max(t^+, t^-) = 2$.

// #remark[
//   The *unsigned* incidence graph $K_(3,3)$ has rank-width 1 (the all-ones cut matrix has rank 1). However, the *signed* rank-width is 2 because $M^+$ and $M^-$ encode different edge patterns. The algorithm's complexity depends on $t = 2$, not the unsigned rank-width.
// ]

// === Step 3: Subspace Structure

// With $t = 2$, the state components are subspaces of $bb(F)_2^2$. The lattice of subspaces is:

// #figure(
//   table(
//     columns: 3,
//     stroke: 0.5pt,
//     table.header([*Dimension*], [*Subspaces*], [*Count*]),
//     [0], [${bold(0)}$], [1],
//     [1], [$chevron.l (1,0) chevron.r, chevron.l (0,1) chevron.r, chevron.l (1,1) chevron.r$], [3],
//     [2], [$bb(F)_2^2$], [1],
//   ),
//   caption: [The 5 subspaces of $bb(F)_2^2$.],
// )

// === Step 4: Labeling

// The labeling assigns vectors in $bb(F)_2^2$ to each vertex. From the cut matrices:
// - $lambda^+(x_1) = (1, 0, 1)^top$, $lambda^-(x_1) = (0, 1, 0)^top$ (column vectors of $M^+, M^-$)
// - $lambda^+(x_2) = lambda^+(x_3) = (1, 1, 0)^top$, $lambda^-(x_2) = lambda^-(x_3) = (0, 0, 1)^top$

// For clauses, the labels come from a basis for the row space. With $t = 2$, we can use:
// - $lambda(C_1) = (1, 0)$, $lambda(C_2) = (0, 1)$ (basis vectors)
// - $lambda(C_3)$ depends on the specific decomposition

// === Step 5: Process Variable Leaves

// At each variable leaf, state $(Sigma^+, Sigma^-, Pi^+, Pi^-)$ tracks which subspaces are generated.

// *Leaf $x_1$*:
// - $x_1 = 1$: generates $Sigma^+ = chevron.l lambda^+(x_1) chevron.r$, $Sigma^- = {bold(0)}$
// - $x_1 = 0$: generates $Sigma^+ = {bold(0)}$, $Sigma^- = chevron.l lambda^-(x_1) chevron.r$

// *Leaves $x_2, x_3$* (identical labels):
// - $x_i = 1$: generates $Sigma^+ = chevron.l lambda^+(x_i) chevron.r$, $Sigma^- = {bold(0)}$
// - $x_i = 0$: generates $Sigma^+ = {bold(0)}$, $Sigma^- = chevron.l lambda^-(x_i) chevron.r$

// === Step 6: Process Clause Leaves

// Each clause $C_j$ requires that its label not be orthogonal to the expectation subspace corresponding to its literals' polarity.

// *Leaf $C_1$* (all positive literals): Requires $lambda(C_1) in.not (Pi^+)^perp$, i.e., $Pi^+$ must contain a vector not orthogonal to $lambda(C_1)$.

// *Leaf $C_2$* (mixed): Has positive literals $x_2, x_3$ and negative literal $x_1$. Satisfaction can come from either $Pi^+$ or $Pi^-$.

// *Leaf $C_3$* (mixed): Has positive literal $x_1$ and negative literals $x_2, x_3$.

// === Step 7: Join and Final Count

// The join operation combines subspaces via sum. At the root, we require $Pi^+ = Pi^- = {bold(0)}$ (all expectations satisfied).

// #figure(
//   table(
//     columns: 4,
//     stroke: 0.5pt,
//     table.header([*Assignment $(x_1, x_2, x_3)$*], [*$Sigma^+$*], [*$Sigma^-$*], [*Satisfies?*]),
//     [$(0, 0, 0)$], [${bold(0)}$], [$bb(F)_2^2$], [No: $C_1$ unsatisfied],
//     [$(1, 0, 0)$], [1-dim], [1-dim], [No: $C_2$ unsatisfied],
//     [$(0, 1, 0)$], [1-dim], [1-dim], [Yes],
//     [$(0, 0, 1)$], [1-dim], [1-dim], [Yes],
//     [$(1, 1, 0)$], [$bb(F)_2^2$], [1-dim], [Yes],
//     [$(1, 0, 1)$], [$bb(F)_2^2$], [1-dim], [Yes],
//     [$(0, 1, 1)$], [1-dim], [$bb(F)_2^2$], [No: $C_3$ unsatisfied],
//     [$(1, 1, 1)$], [$bb(F)_2^2$], [${bold(0)}$], [Yes],
//   ),
//   caption: [Final satisfiability via subspace coverage. "1-dim" denotes a 1-dimensional subspace.],
// )

// The expectation matching verifies:
// - $C_1$ (all positive): needs $Sigma^+$ to cover $lambda(C_1)$
// - $C_2$ (positive $x_2, x_3$; negative $x_1$): needs $Sigma^+$ or $Sigma^-$ to cover $lambda(C_2)$
// - $C_3$ (positive $x_1$; negative $x_2, x_3$): needs $Sigma^+$ or $Sigma^-$ to cover $lambda(C_3)$

// $ \#op("SAT")(phi) = 5 $

// #remark[
//   With signed rank-width $t = 2$, the DP tracks states over 5 subspaces per component, yielding $O(5^4) = 625$ possible states per node. This is larger than the $2^3 = 8$ assignments but demonstrates the general structure. For formulas where $t << n$, the compression becomes significant.
// ]

// == Complexity Analysis

// #figure(
//   table(
//     columns: 3,
//     stroke: 0.5pt,
//     table.header([*Component*], [*Size*], [*For rank-width $t$*]),
//     [Subspaces of $bb(F)_2^t$], [$S(t) <= 2^(t(t+1)\/4)$], [Lemma 3.3 @ganian2010],
//     [Table size per node], [$O(S(t)^4)$], [4 subspaces per state],
//     [Join operation], [$O(t^3 dot S(t)^6)$], [Loop over 6 subspaces],
//   ),
//   caption: [DP complexity depends on the signed rank-width $t$.],
// ) <tab:dp-complexity>

// The actual complexity from @ganian2010 is:
// $ O(t^3 dot 2^(3t(t+1)\/2) dot |phi|) $
// where $t = max(t^+, t^-)$ is the signed rank-width. Although there are $O(S(t)^4)$ possible states, the join operation can be computed efficiently by iterating over 6 subspaces (the remaining 6 are determined by linear algebra), with each iteration requiring $O(t^3)$ time.


// = Computing Rank-Width

// Computing rank-width exactly is NP-hard @oum2006. Several algorithmic approaches exist:

// #v(10pt)
// #figure(
//   table(
//     columns: 3,
//     stroke: 0.5pt,
//     table.header([*Approach*], [*Complexity*], [*Notes*]),
//     [Exact], [$O(3^n)$], [Brute force],
//     [FPT], [$O(f(k) dot n^3)$], [Fixed-parameter tractable],
//     [Approximation], [Polynomial], [$(3k+1)$-approximation @oum2006],
//   ),
//   caption: [Algorithms for computing rank-width.],
// )

// The cut-rank function $rho(A) = op("rank")_(bb(F)_2)(M_A)$ satisfies two key properties:
// - *Symmetric*: $rho(A) = rho(B)$
// - *Submodular*: $rho(A) + rho(B) >= rho(A inter B) + rho(A union B)$

// #figure(
//   rect(
//     width: 100%,
//     inset: 12pt,
//     fill: rgb("#f0fff0"),
//     stroke: rgb("#4a4"),
//     [
//       *Why submodularity helps*: Submodular functions have "diminishing returns" --- adding an element to a larger set increases the function value by at most as much as adding it to a smaller set. This property enables:

//       1. *Greedy algorithms*: Local choices lead to globally good solutions
//       2. *Branch decomposition*: The $(3k+1)$-approximation @oum2006 uses submodularity to efficiently explore the space of decompositions
//       3. *Polynomial-time verification*: Given a decomposition, we can verify its width in polynomial time

//       Without submodularity, finding optimal decompositions would require exhaustive search over all possible tree structures.
//     ]
//   ),
// )

// = Complexity Comparison

// #figure(
//   table(
//     columns: 3,
//     stroke: 0.5pt,
//     table.header([*Algorithm*], [*Time Complexity*], [*Reference*]),
//     [Treewidth-based], [$O(2^(t w) dot n)$], [Standard DP],
//     [Clique-width-based], [$O(2^(c w) dot n^2)$], [@courcelle2000],
//     [Rank-width-based], [$O(2^(3 r w) dot n^2 dot m)$], [@ganian2010],
//   ),
//   caption: [Parameterized algorithms for \#SAT.],
// ) <tab:complexity>


= Background Notes

The algorithmic exploitation of neighborhood equivalence over finite fields emerged from two parallel developments in parameterized complexity.

*The $d$-neighbor equivalence framework.* Bui-Xuan, Telle, and Vatshelle @bui-xuan2013 introduced the $d$-neighbor equivalence relation for dynamic programming on graph decompositions. Their key insight was that for locally checkable problems (where feasibility depends only on local neighborhoods), it suffices to track equivalence classes rather than individual vertex subsets. This yields polynomial-time algorithms for graphs of bounded boolean-width @bui-xuan2011, with runtime $O(n^4 2^(O(k^2)))$ for decompositions of boolean-width $k$.

*The rank-based approach.* Independently, Bodlaender, Cygan, Kratsch, and Nederlof @bodlaender2015 developed the rank-based approach for connectivity problems parameterized by treewidth. Their central observation was that *representative sets*---small subsets that preserve all relevant information for dynamic programming---can be computed efficiently via Gaussian elimination over finite fields. This technique achieves single-exponential runtime $2^(O(k)) n^(O(1))$ for problems like Steiner Tree and Hamiltonian Path, where previous algorithms required $k^(O(k)) n^(O(1))$ time.

*Synthesis for rank-width.* The Ganian--Hliněný--Obdržálek algorithm @ganian2010 synthesizes these ideas for \#SAT parameterized by rank-width. The cut-rank function $rho_G (X) = op("rank")_(bb(F)_2)(M_X)$ directly measures the number of equivalence classes, and the use of $bb(F)_2$ enables XOR convolution for efficient joins. The signed graph formulation extends this to track literal polarity, with the expectation mechanism (originating from @bui-xuan2011 for dominating set) enabling single-exponential dependence on the signed rank-width.

#bibliography("refs.bib", title: "References", style: "ieee")
