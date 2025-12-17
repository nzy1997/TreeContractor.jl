#import "@preview/clear-iclr:0.7.0": *
#import "@preview/cetz:0.4.2": canvas, draw

#show: iclr.with(
  title: [Rank-Width for \#SAT],
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
    We provide a pedagogical introduction to rank-width and its application to the propositional model counting problem \#SAT. Through a concrete example, we demonstrate how rank-width can lead to exponential speedups over treewidth-based approaches for certain graph families.
  ],
  accepted: none,
  bibliography: none,
)

= Introduction

Ganian, Hliněný, and Obdržálek @ganian2010 provide a \#SAT algorithm with runtime single-exponential in the _rank-width_ of the incidence graph. The key relationship is $op("rw")(G) <= op("tw")(G) + 1$, but rank-width can be exponentially smaller for certain graph families.

= Rank Decomposition and Rank-Width

A *rank decomposition* $(T, L)$ of graph $G = (V, E)$ consists of a subcubic tree $T$ and a bijection $L$ from leaves to $V$. Each edge $e in T$ induces a partition $(A_e, B_e)$ of $V$ ($B_e := V - A_e$).

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
  caption: [A rank decomposition. Removing edge $e$ partitions the leaves into $A_e$ and $B_e$.],
) <fig:rank-decomp>

#v(10pt)
For partition $(A, B)$, the *cut matrix* $M_A in bb(F)_2^(|A| times |B|)$ has $M_A [u,v] = 1$ iff $(u,v) in E$. The *rank-width* is:
$ op("rw")(G) = min_((T,L)) max_(e in T) op("rank")_(bb(F)_2)(M_(A_e)) $

= Neighborhood Equivalence

Let $G = (V, E)$ be a graph. For any partition $(A, B)$ of $V$ and any $Y subset.eq A$, define:
$ N_B (Y) = {v in B : exists y in Y, {y, v} in E} $
the set of vertices in $B$ adjacent to at least one vertex in $Y$.

Two subsets $Y, Y' subset.eq A$ are *neighborhood-equivalent* w.r.t. the partition $(A, B)$ if $N_B (Y) = N_B (Y')$.

Note: The equivalence relation depends on the choice of partition. Each edge $e$ of the rank decomposition induces a partition $(A_e, B_e)$, and thus a different equivalence relation on subsets of $A_e$.

The number of equivalence classes for partition $(A, B)$ is at most $2^(op("rank")(M_A))$, since each class corresponds to a distinct vector in the row space of $M_A$ @bui-xuan2011. This enables DP compression: instead of $2^(|A|)$ subsets, we track at most $2^(op("rw")(G))$ classes at each node.

#figure(
  canvas(length: 1cm, {
    import draw: *
    let inside-col = rgb("#4a90d9")
    let boundary-col = rgb("#e07b53")
    let class1-col = rgb("#9b59b6")  // purple for class 1
    let class2-col = rgb("#27ae60")  // green for class 2
    // Inside region
    rect((-4, -2.5), (1, 2.5), fill: inside-col.lighten(90%), stroke: inside-col)
    content((-1.5, 2.8), text(9pt, inside-col)[$A_e$])
    // Boundary region
    rect((1.5, -2.5), (4.5, 2.5), fill: boundary-col.lighten(90%), stroke: boundary-col)
    content((3, 2.8), text(9pt, boundary-col)[$B_e$])
    // Vertices inside - 4 vertices
    circle((-3, 1.5), radius: 0.3, fill: class1-col.lighten(50%), stroke: class1-col + 1.5pt, name: "a1")
    content((-3, 1.5), text(8pt)[$a_1$])
    circle((0, 2.8), radius: 0.3, fill: class1-col.lighten(50%), stroke: class1-col + 1.5pt, name: "a2")
    content((0, 2.8), text(8pt)[$a_2$])
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
  caption: [Four vertices in $A_e$, but only 2 equivalence classes based on their neighborhood in $B_e$.],
) <fig:cut-boundary>

In @fig:cut-boundary, vertices $a_1, a_2, a_3$ (purple) all connect only to $b_1$, so they share the same neighborhood $N_B = {b_1}$ and belong to *Class 1*. Vertex $a_4$ (green) connects to both $b_1$ and $b_2$, so $N_B ({a_4}) = {b_1, b_2}$ --- *Class 2*.

With 4 vertices, there are $2^4 = 16$ possible subsets of $A_e$. But since $a_1, a_2, a_3$ are interchangeable (same neighborhood), many subsets collapse:
- ${a_1}, {a_2}, {a_3}$ all have neighborhood ${b_1}$ --- same class
- ${a_1, a_2}, {a_1, a_3}, {a_2, a_3}$ all have neighborhood ${b_1}$ --- same class

== Example: Dense Formula --- Where Neighborhood Equivalence Helps

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

Recall from the neighborhood equivalence definition that two sets $S, S' subset.eq X$ are equivalent iff $N(S) inter B = N(S') inter B$. To compute this efficiently, we represent the neighborhood as a *characteristic vector* over $bb(F)_2$.

#figure(
  rect(
    width: 100%,
    inset: 12pt,
    fill: rgb("#f0f8ff"),
    stroke: rgb("#4a90d9"),
    [
      *Definition (Neighborhood Characteristic)*: Given a set $S subset.eq X$ with characteristic vector $bold(v) in bb(F)_2^n$ (where $v_i = 1$ iff $x_i in S$), the *neighborhood characteristic* is:
      $ bold(n)_S = bold(v)^top M in bb(F)_2^m $
      where $M$ is the cut matrix. The $j$-th entry is:
      $ (bold(n)_S)_j = plus.big_(i : x_i in S) M[i,j] mod 2 $
      This is the characteristic vector of $N(S) inter B$ computed over $bb(F)_2$.
    ]
  ),
  caption: [The neighborhood characteristic is the $bb(F)_2$-representation of the neighborhood.],
)
#v(10pt)
*For \#SAT*: When $S$ is the set of TRUE variables and $B$ is the set of clauses, $(bold(n)_S)_j$ counts (mod 2) how many TRUE variables appear in clause $C_j$.

*Connection to neighborhood equivalence*: Two sets $S, S'$ are neighborhood-equivalent iff they have the same neighborhood characteristic:
$ S tilde.eq S' <==> bold(n)_S = bold(n)_(S') <==> bold(v)^top M = bold(v')^top M $

=== Computing the Neighborhood Characteristic

Let $bold(v) in bb(F)_2^3$ be the characteristic vector of TRUE variables (e.g., $bold(v) = (1,1,0)$ means $x_1 = x_2 = "TRUE", x_3 = "FALSE"$).

For our dense formula with $M = J$ (all-ones matrix):
$ bold(n)_S = bold(v)^top M = bold(v)^top mat(1, 1, 1; 1, 1, 1; 1, 1, 1) = (v_1 xor v_2 xor v_3) dot (1, 1, 1) $

Since every variable appears in every clause, the neighborhood characteristic depends only on the *parity* of $|S|}$!

#v(10pt)
#figure(
  table(
    columns: 5,
    stroke: 0.5pt,
    table.header([*TRUE vars $S$*], [*$bold(v)$*], [*$|S| mod 2$*], [*Neighborhood $bold(n)_S$*], [*Class*]),
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

=== Why is $\{x_1, x_2\}$ in Class 0?

Let's trace through the computation step by step for assignment $S = {x_1, x_2}$ (i.e., $x_1 = x_2 = "TRUE", x_3 = "FALSE"$):

1. *Characteristic vector*: $bold(v) = (1, 1, 0)$

2. *Matrix multiplication*:
$ bold(n)_S = bold(v)^top M = mat(1, 1, 0) mat(1, 1, 1; 1, 1, 1; 1, 1, 1) $

3. *For each clause $C_j$*, compute $(bold(n)_S)_j = v_1 dot M[1,j] xor v_2 dot M[2,j] xor v_3 dot M[3,j]$:
   - $(bold(n)_S)_1 = 1 dot 1 xor 1 dot 1 xor 0 dot 1 = 1 xor 1 xor 0 = 0$
   - $(bold(n)_S)_2 = 1 dot 1 xor 1 dot 1 xor 0 dot 1 = 1 xor 1 xor 0 = 0$
   - $(bold(n)_S)_3 = 1 dot 1 xor 1 dot 1 xor 0 dot 1 = 1 xor 1 xor 0 = 0$

4. *Result*: $bold(n)_S = (0, 0, 0)$ --- Class 0!

Since both $x_1$ and $x_2$ appear in every clause, when both are TRUE, each clause contributes $1 xor 1 = 0$ (even parity). This yields the same neighborhood characteristic as $emptyset$.

=== Implications for \#SAT

Two assignments are *neighborhood-equivalent* if they produce the same neighborhood characteristic. Equivalent assignments interact identically with the boundary vertices during the DP computation.

#figure(
  rect(
    width: 100%,
    inset: 12pt,
    fill: rgb("#f8fff8"),
    stroke: rgb("#4a4"),
    [
      *Implication*: When counting satisfying assignments, we only need to track *how many* assignments fall into each equivalence class, not which specific assignments. This reduces the state space from $2^n$ to $2^k$ where $k$ is the rank-width.
    ]
  ),
  caption: [Equivalence classes enable efficient counting.],
)

#v(10pt)
In this example, only 2 distinct neighborhood characteristics exist:
- *Class 0* (neighborhood $(0,0,0)$): $emptyset, {x_1,x_2}, {x_1,x_3}, {x_2,x_3}$
- *Class 1* (neighborhood $(1,1,1)$): ${x_1}, {x_2}, {x_3}, {x_1,x_2,x_3}$

This yields a 4x compression: 8 assignments reduce to 2 equivalence classes.

== Counterexample: Sparse Formula --- No Compression

Now consider a CNF formula where *each variable appears in exactly one clause*:
$ psi = C_1 and C_2 and C_3 $
where:
- $C_1 = (x_1)$
- $C_2 = (x_2)$
- $C_3 = (x_3)$

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
The cut matrix separating variables from clauses is the *identity matrix*:
$ M = mat(1, 0, 0; 0, 1, 0; 0, 0, 1) $

Over $bb(F)_2$, this has *rank 3* (full rank). Let's compute the neighborhood characteristics:

For any set $S$ with characteristic vector $bold(v)$:
$ bold(n)_S = bold(v)^top I = bold(v) $

The neighborhood characteristic *equals* the characteristic vector itself! This means every distinct set has a distinct neighborhood.

#v(10pt)
#figure(
  table(
    columns: 4,
    stroke: 0.5pt,
    table.header([*TRUE variables $S$*], [*$bold(v)$*], [*Neighborhood $bold(n)_S = bold(v)^top I$*], [*Class*]),
    [$emptyset$], [$(0,0,0)$], [$(0,0,0)$], [$(0,0,0)$],
    [${x_1}$], [$(1,0,0)$], [$(1,0,0)$], [$(1,0,0)$],
    [${x_2}$], [$(0,1,0)$], [$(0,1,0)$], [$(0,1,0)$],
    [${x_3}$], [$(0,0,1)$], [$(0,0,1)$], [$(0,0,1)$],
    [${x_1, x_2}$], [$(1,1,0)$], [$(1,1,0)$], [$(1,1,0)$],
    [${x_1, x_3}$], [$(1,0,1)$], [$(1,0,1)$], [$(1,0,1)$],
    [${x_2, x_3}$], [$(0,1,1)$], [$(0,1,1)$], [$(0,1,1)$],
    [${x_1, x_2, x_3}$], [$(1,1,1)$], [$(1,1,1)$], [$(1,1,1)$],
  ),
  caption: [Each truth assignment has a *unique* neighborhood characteristic --- no compression.],
) <tab:matching-equiv>

*Key observation*: Every truth assignment has a *distinct* neighborhood characteristic! The 8 assignments map to 8 different equivalence classes.

*Why?* Each variable $x_i$ only appears in clause $C_i$. Setting $x_i$ to TRUE affects *only* $C_i$, independently of other variables. The clauses can distinguish every assignment because no two variables "overlap" in their clause memberships.

#figure(
  rect(
    width: 100%,
    inset: 12pt,
    fill: rgb("#fff0f0"),
    stroke: rgb("#c44"),
    [
      *Why sparse formulas don't benefit*: When each variable appears in a unique clause, the cut matrix is the identity (full rank). Every truth assignment produces a unique neighborhood characteristic:
      $ "equivalence classes" = 2^(op("rank")(I_n)) = 2^n = "number of assignments" $
      No compression from neighborhood equivalence!
    ]
  ),
  caption: [Sparse variable-clause structure prevents algebraic compression.],
)

== Comparison: Why Formula Structure Matters

#figure(
  table(
    columns: 5,
    stroke: 0.5pt,
    table.header([*Formula Type*], [*Cut Matrix*], [*Rank*], [*Equiv. Classes*], [*Compression*]),
    [Dense (every var in every clause)], [All-ones $J$], [$1$], [$2^1 = 2$], [$2^(n-1) times$],
    [Sparse (each var in one clause)], [Identity $I$], [$n$], [$2^n$], [None],
  ),
  caption: [Dense formulas enable compression; sparse formulas do not.],
) <tab:compression-comparison>

The contrast is clear: dense formulas (where variables share clauses) have low-rank cut matrices, while sparse formulas (where variables appear in disjoint clauses) have full-rank cut matrices. This determines whether neighborhood equivalence provides any benefit for \#SAT.

== Scaling to $K_(n,n)$

The $K_(3,3)$ example generalizes to $K_(n,n)$:
#v(10pt)
#figure(
  table(
    columns: 4,
    stroke: 0.5pt,
    table.header([*Graph*], [*Treewidth*], [*Rank-width*], [*Compression*]),
    [$K_(n,n)$], [$n$], [$1$], [$2^(n-1) times$],
  ),
  caption: [Complete bipartite graphs have constant rank-width but linear treewidth.],
) <tab:knn>

For $n = 100$: treewidth-based algorithms require $2^(100) approx 10^(30)$ states, while rank-width-based algorithms require only $2^1 = 2$ states.

= The DP Algorithm

The Ganian et al. algorithm @ganian2010 processes the rank decomposition tree bottom-up, using neighborhood equivalence to compress the state space.

== Overview

The algorithm traverses the rank decomposition $(T, L)$ from leaves to root. At each node, it maintains a table mapping *neighborhood classes* to *model counts*.

#figure(
  rect(
    width: 100%,
    inset: 12pt,
    fill: rgb("#fff8e8"),
    stroke: rgb("#c90"),
    [
      *Key Idea*: For a subtree with vertex set $A$ and boundary $B := V - A$:
      - Group partial assignments by their neighborhood equivalence class
      - Track how many satisfying assignments fall into each class
      - At most $2^(rho_G (A))$ classes to track (not $2^(|A|)$)
    ]
  ),
  caption: [The compression principle behind the DP algorithm.],
)

== State Representation

At each edge $e$ of the decomposition tree (separating $A_e$ from $B_e$):
- *State*: A vector $bold(s) in bb(F)_2^(|B_e|)$ representing the neighborhood characteristic on the boundary
- *Table*: $T_e [bold(s)]$ = number of partial truth assignments to variables in $A_e$ that:
  1. Satisfy all clauses fully contained in $A_e$
  2. Have neighborhood characteristic $bold(s)$ on the boundary

Since equivalent neighborhood characteristics collapse, the table has at most $2^(rho_G (A_e))$ entries.

== Processing Nodes

The algorithm processes three types of nodes:

=== Leaf Nodes (Variable Introduction)

For a leaf corresponding to variable $v$ with boundary $B$:
- Compute $N(v) inter B$ (neighbors of $v$ in the boundary)
- $T[bold(0)]$ += 1 (assignment $v = 0$, no neighborhood contribution)
- $T[chi_(N(v) inter B)]$ += 1 (assignment $v = 1$, contributes its neighborhood)

where $chi_S$ is the characteristic vector of set $S$.

=== Internal Nodes (Join)

For an internal node with children corresponding to edges $e_1, e_2$:
$ T_e [bold(s)] = sum_(bold(s)_1 xor bold(s)_2 = bold(s)) T_(e_1)[bold(s)_1] dot T_(e_2)[bold(s)_2] $

This is a *convolution over $bb(F)_2$* --- the XOR reflects that neighborhoods combine additively in $bb(F)_2$.

=== Root Node (Final Count)

At the root, sum over all states that correspond to *satisfied* formulas:
$ \#op("SAT")(phi) = sum_(bold(s) : "all clauses satisfied") T_("root")[bold(s)] $

== Complexity Analysis

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([*Component*], [*Size*], [*For rank-width $k$*]),
    [Equivalence classes], [$2^(rho_G (A_e))$], [$<= 2^k$],
    [Table size per node], [$O(2^k)$], [$O(2^k)$],
    [Join operation], [$O(|T|^2)$], [$O(2^(2k))$],
    [Total (naive)], [$O(n dot 2^(2k))$], [$O(n dot 2^(2k))$],
  ),
  caption: [DP complexity depends on the rank-width $k$.],
) <tab:dp-complexity>

The actual complexity from @ganian2010 is $O(2^(3k) dot n^2 dot m)$ where $n$ is the number of variables and $m$ is the number of clauses.

= Additional Graph Families

#figure(
  table(
    columns: 4,
    stroke: 0.5pt,
    table.header([*Graph Family*], [*Treewidth*], [*Rank-width*], [*Advantage*]),
    [$K_(n,n)$ (complete bipartite)], [$n$], [$1$], [Exponential],
    [$K_n$ (complete)], [$n-1$], [$ceil(n/2) - 1$], [None],
    [$P_n$ (path)], [$1$], [$1$], [None],
    [$C_n$ (cycle)], [$2$], [$2$], [None],
    [Grid $G_(m times n)$], [$min(m,n)$], [$Theta(min(m,n))$], [None],
  ),
  caption: [Comparison of treewidth and rank-width for various graph families.],
) <tab:graph-comparison>

Rank-width captures algebraic structure (low-rank adjacency matrices) rather than combinatorial sparsity. Dense graphs with regular structure can have low rank-width.

= Computing Rank-Width

Computing rank-width exactly is NP-hard @oum2006. Several algorithmic approaches exist:

#v(10pt)
#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([*Approach*], [*Complexity*], [*Notes*]),
    [Exact], [$O(3^n)$], [Brute force],
    [FPT], [$O(f(k) dot n^3)$], [Fixed-parameter tractable],
    [Approximation], [Polynomial], [$(3k+1)$-approximation @oum2006],
  ),
  caption: [Algorithms for computing rank-width.],
)

The cut-rank function $rho(A) = op("rank")_(bb(F)_2)(M_A)$ satisfies two key properties:
- *Symmetric*: $rho(A) = rho(B)$
- *Submodular*: $rho(A) + rho(B) >= rho(A inter B) + rho(A union B)$

Submodularity enables efficient approximation algorithms.

= Complexity Comparison

#figure(
  table(
    columns: 3,
    stroke: 0.5pt,
    table.header([*Algorithm*], [*Time Complexity*], [*Reference*]),
    [Treewidth-based], [$O(2^(t w) dot n)$], [Standard DP],
    [Clique-width-based], [$O(2^(c w) dot n^2)$], [@courcelle2000],
    [Rank-width-based], [$O(2^(3 r w) dot n^2 dot m)$], [@ganian2010],
  ),
  caption: [Parameterized algorithms for \#SAT.],
) <tab:complexity>

= Conclusion

Rank-width provides an alternative parameterization for \#SAT. By measuring the $bb(F)_2$-rank of cut matrices rather than separator size, it can yield exponential improvements for graph families with low-rank algebraic structure.

#bibliography("refs.bib", title: "References", style: "ieee")
