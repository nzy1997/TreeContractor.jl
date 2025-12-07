#import "@preview/cetz:0.4.0": canvas, draw, tree, coordinate
#import "@preview/cetz-plot:0.1.2": *
#import "@preview/ctheorems:1.1.3": *

#set math.equation(numbering: "(1)")
#show link: set text(blue)
#show heading.where(level: 1): set text(20pt)
#show: thmrules

#show raw.where(block: true): it=>{
  block(fill:rgb("#fcf9ec"),inset:1.5em,width:99%,text(it))
}

#let labelnode(loc, label, name: none) = {
  import draw: *
  content(loc, text(black, label), align: center, fill:silver, frame:"rect", padding:0.07, stroke: none, name: name)
}
#let labeledge(from, to, label, name: none) = {
  import draw: *
  line(from, to, name: "line")
  if label != none {
    labelnode("line.mid", label, name: name)
  }
}

#let tensor(location, name, label) = {
  import draw: *
  circle(location, radius: 10pt, name: name)
  content((), text(black, label))
}

#let definition = thmbox("definition", "Definition", inset: (x: 1.2em, top: 1em, bottom: 1em), base: none, stroke: none, fill: rgb("#e8f4fd"), namefmt: x => [(#strong[#x.])], titlefmt: x => [(#emph[#x])])
#let theorem = thmbox("theorem", "Theorem", base: none, stroke: none, fill: rgb("#f0f9e8"), namefmt: x => [(#strong[#x.])], titlefmt: x => [(#emph[#x])])
#let lemma = thmbox("lemma", "Lemma", base: "theorem", stroke: none, fill: rgb("#f0f9e8"), namefmt: x => [(#strong[Lemma #x.])], titlefmt: x => [(#emph[#x])])
#let corollary = thmbox("corollary", "Corollary", base: "theorem", stroke: none, fill: rgb("#f0f9e8"), namefmt: x => [(#strong[Corollary #x.])], titlefmt: x => [(#emph[#x])])
#let proposition = thmbox("proposition", "Proposition", base: "theorem", stroke: none, fill: rgb("#f0f9e8"), namefmt: x => [(#strong[Proposition #x.])], titlefmt: x => [(#emph[#x])])
#let proof = thmproof("proof", "Proof")
#let ket(it) = [$|#it angle.r$]


= Tree contractor

Given a tree decomposition of a tensor network, can we contract the tensor network efficiently?

+ Determine a tree decomposition of the tensor network (@fig:tree-decomposition (c)).
  #figure(canvas({
  import draw: *
  let d = 1.1
  let s(it) = text(11pt, it)
  let locs_labels = ((0, 0), (d, 0), (0, -d), (0, -2 * d), (d, -2 * d), (2 * d, 0), (2 * d, -d), (2 * d, -2 * d))
  for (loc, t, name) in (((0.5 * d, -0.5 * d), s[$T_1$], "T_1"), ((1.5 * d, -0.5 * d), s[$T_2$], "T_2"), ((1.5 * d, -1.5 * d), s[$T_3$], "T_3"), ((0.5 * d, -1.5 * d), s[$T_4$], "T_4")) {
    circle(loc, radius: 0.3, name: name)
    content(loc, s[#t])
  }
  for ((loc, t), name) in locs_labels.zip((s[$A$], s[$B$], s[$C$], s[$D$], s[$E$], s[$F$], s[$G$], s[$H$])).zip(("A", "B", "C", "D", "E", "F", "G", "H")) {
    labelnode(loc, t, name: name)
  }
  for (src, dst) in (("A", "T_1"), ("B", "T_1"), ("C", "T_1"), ("F", "T_2"), ("G", "T_2"), ("B", "T_2"), ("H", "T_3"), ("E", "T_3"), ("G", "T_3"), ("D", "T_4"), ("C", "T_4"), ("E", "T_4")) {
    line(src, dst)
  }
  content((d, -3), text(12pt)[(a)])
  content((3.5, -1), text(12pt)[$arrow.double.r$])
  content((3.5, -1.5), text(10pt)[Line graph])
  set-origin((5, 0))
  let colors = (color.hsv(30deg, 90%, 70%), color.hsv(120deg, 90%, 70%), color.hsv(210deg, 90%, 70%), color.hsv(240deg, 90%, 70%), color.hsv(330deg, 90%, 70%), color.hsv(120deg, 90%, 70%), color.hsv(210deg, 90%, 70%), color.hsv(240deg, 90%, 70%))
  let texts = ("A", "B", "C", "D", "E", "F", "G", "H")
  for (loc, color, t) in locs_labels.zip(colors, texts) {
    circle(loc, radius: 0.3, name: t)
    content(loc, text(12pt, color)[#t])
  }
  for (a, b) in (("A", "B"), ("A", "C"), ("B", "C"), ("C", "D"), ("C", "E"), ("D", "E"), ("E", "G"), ("G", "H"), ("E", "H"), ("F", "G"), ("F", "B"), ("B", "G")) {
    line(a, b)
  }
  content((d, -3), text(12pt)[(b)])
  content((3.5, -1), text(12pt)[$arrow.double.r$])
  content((3.5, -1.5), text(10pt)[T. D.])
  set-origin((5, 0))
  for (loc, bag) in (((0, 0), "B1"), ((0, -2), "B2"), ((1, -1), "B3"), ((3, -1), "B4"), ((4, 0), "B5"), ((4, -2), "B6")) {
    circle(loc, radius: 0.55, name: bag)
    content((rel: (0, -0.75)), text(10pt, gray)[#bag])
  }
  let topleft = (-0.2, 0.2)
  let topright = (0.2, 0.2)
  let bottom = (0, -0.3)
  let top = (0, 0.3)
  let bottomleft = (-0.2, -0.2)
  let bottomright = (0.2, -0.2)
  let right = (0.3, 0)
  let left = (-0.3, 0)
  content((rel:topright, to: "B1"), text(10pt, colors.at(1))[B], name: "b1")
  content((rel:topleft, to: "B1"), text(10pt, colors.at(0))[A], name: "a1")
  content((rel:bottom, to: "B1"), text(10pt, colors.at(2))[C], name: "c1")

  content((rel:top, to: "B2"), text(10pt, colors.at(2))[C], name: "c2")
  content((rel:bottomleft, to: "B2"), text(10pt, colors.at(3))[D], name: "d1")
  content((rel:right, to: "B2"), text(10pt, colors.at(4))[E], name: "e1")

  content((rel:topright, to: "B3"), text(10pt, colors.at(1))[B], name: "b2")
  content((rel:left, to: "B3"), text(10pt, colors.at(2))[C], name: "c3")
  content((rel:bottomright, to: "B3"), text(10pt, colors.at(4))[E], name: "e2")

  content((rel:topleft, to: "B4"), text(10pt, colors.at(1))[B], name: "b3")
  content((rel:bottomleft, to: "B4"), text(10pt, colors.at(4))[E], name: "e3")
  content((rel:right, to: "B4"), text(10pt, colors.at(6))[G], name: "g1")

  content((rel:left, to: "B5"), text(10pt, colors.at(1))[B], name: "b4")
  content((rel:topright, to: "B5"), text(10pt, colors.at(5))[F], name: "f1")
  content((rel:bottom, to: "B5"), text(10pt, colors.at(6))[G], name: "g2")

  content((rel:left, to: "B6"), text(10pt, colors.at(4))[E], name: "e4")
  content((rel:top, to: "B6"), text(10pt, colors.at(6))[G], name: "g3")
  content((rel:bottomright, to: "B6"), text(10pt, colors.at(7))[H], name: "h1")

  line("b1", "b2", stroke: colors.at(1))
  line("b2", "b3", stroke: colors.at(1))
  line("b3", "b4", stroke: colors.at(1))
  line("c1", "c3", stroke: colors.at(2))
  line("c2", "c3", stroke: colors.at(2))
  line("e1", "e2", stroke: colors.at(4))
  line("e2", "e3", stroke: colors.at(4))
  line("e3", "e4", stroke: colors.at(4))
  line("g1", "g2", stroke: colors.at(6))
  line("g1", "g3", stroke: colors.at(6))
  content((2, -3), text(12pt)[(c)])
}),
caption: [(a) A tensor network. (b) A line graph for the tensor network. Labels are connected if and only if they appear in the same tensor. (c) A tree decomposition (T. D.) of the line graph.]
) <fig:tree-decomposition>


+ Assign each label a arbitrary tree node that contains the label, each tensor to a tree bag that contains all the labels of the tensor.

  #figure(canvas({
  import draw: *
  let colors = (color.hsv(30deg, 90%, 70%), color.hsv(120deg, 90%, 70%), color.hsv(210deg, 90%, 70%), color.hsv(240deg, 90%, 70%), color.hsv(330deg, 90%, 70%), color.hsv(120deg, 90%, 70%), color.hsv(210deg, 90%, 70%), color.hsv(240deg, 90%, 70%))
  for (loc, name, bag) in (((0, 0), "T1", [$T_1$]), ((0, -2), "T4", [$T_4$]), ((1, -1), "M1", []), ((3, -1), "M2", []), ((4, 0), "T2", [$T_2$]), ((4, -2), "T3", [$T_3$])) {
    circle(loc, radius: 0.55, name: name)
    content((rel: (0, -0.75)), text(10pt, gray)[#bag])
  }
  let topleft = (-0.2, 0.2)
  let topright = (0.2, 0.2)
  let bottom = (0, -0.3)
  let top = (0, 0.3)
  let bottomleft = (-0.2, -0.2)
  let bottomright = (0.2, -0.2)
  let right = (0.3, 0)
  let left = (-0.3, 0)
  content((rel:topleft, to: "T1"), text(10pt, colors.at(0))[A], name: "a")
  content((rel:bottom, to: "T1"), text(10pt, colors.at(2))[C], name: "c")

  content((rel:bottomleft, to: "T4"), text(10pt, colors.at(3))[D], name: "d")

  content((rel:right, to: "M2"), text(10pt, colors.at(6))[G], name: "g")

  content((rel:left, to: "T2"), text(10pt, colors.at(1))[B], name: "b")
  content((rel:topright, to: "T2"), text(10pt, colors.at(5))[F], name: "f")

  content((rel:left, to: "T3"), text(10pt, colors.at(4))[E], name: "e")
  content((rel:bottomright, to: "T3"), text(10pt, colors.at(7))[H], name: "h")

  line("a", "c", stroke: colors.at(2))
  line("c", "M1.center", stroke: colors.at(2))
  line("d", "M1.center", stroke: colors.at(2))
  line("g", "M1.center", stroke: colors.at(2))
  line("a", "c", stroke: colors.at(2))
  line("a", "c", stroke: colors.at(2))
  line("g", "e", stroke: colors.at(2))
  line("e", "h", stroke: colors.at(2))
  line("g", "b", stroke: colors.at(2))
  line("b", "f", stroke: colors.at(2))

  let dy = 0.4
  let cup = (rel: (0, dy), to: "c")
  let dup = (rel: (0, dy), to: "d")
  let eup = (rel: (0, dy), to: "e")
  let gup = (rel: (0, dy), to: "g")
  let mup = (rel: (0, dy), to: "M1")
  line(cup, mup, dup, stroke: red, name: "l1")
  line(eup, gup, mup, dup, stroke: red, name: "l2")
  line(cup, (anchor: 30%, name: "l1"), stroke: red, mark: (end: "straight"))
  line(eup, (anchor: 12%, name: "l2"), stroke: red, mark: (end: "straight"))
  line(gup, (anchor: 45%, name: "l2"), stroke: red, mark: (end: "straight"))
  line(mup, (anchor: 80%, name: "l1"), stroke: red, mark: (end: "straight"))
  line(mup, (anchor: 70%, name: "l1"), stroke: red, mark: (end: "straight"))
  line(cup, "c", stroke: red)
  line(eup, "e", stroke: red)
  line(gup, "g", stroke: red)
  line(dup, "d", stroke: red)
  for (node, color) in ((dup, black), (eup, red), (gup, red), (mup, red), (cup, red)) {
    circle(node, radius: 0.1, fill: color, stroke: none)
  }
}))
+ Perform _tree factorization_ to each tensor. The _tree factorization_ is a tensor factorization determined by the underlying tree structure. e.g. the above $T_4$ can be factored as the MPS/MPO shown in red, where red tensors denotes $delta$ tensors and the black tensors denotes the original tensor $T_4$.

#theorem([The bond dimension of the above tree tensor network does not exceed the tree decomposition's maximum separator size.])

== \#SAT and rank width

The \#SAT (model counting) problem asks: given a Boolean formula in CNF (conjunctive normal form), how many satisfying assignments exist?

For example, the formula $(x_1 or x_2) and (x_2 or x_3) and (x_3 or x_4)$ has multiple satisfying assignments. Model counting is \#P-complete, but becomes tractable for formulas with bounded structural width.

=== Rank width

The _rank width_ is a structural parameter that characterizes the difficulty of solving \#SAT using tree-based decomposition methods. It is defined based on the rank decomposition of the _incidence graph_.

#definition([Incidence graph])[
  For a CNF formula $phi$ with variables $V$ and clauses $C$, the _incidence graph_ $G = (V union C, E)$ is a bipartite graph where each variable $v in V$ is connected to each clause $c in C$ that contains $v$ or $overline(v)$.
]

#definition([Rank width])[
  Let $G = (V union C, E)$ be the incidence graph of a CNF formula. For a bipartition $(A, B)$ of vertices, the _cut-rank_ is:
  $
    "cut-rank"(A, B) = "rank"(M_(A, B))
  $
  where $M_(A, B)$ is the submatrix of the adjacency matrix with rows indexed by $A$ and columns indexed by $B$, and the rank is computed over $bb(F)_2$ (the binary field).
  
  A _rank decomposition_ of $G$ is a binary tree $T$ whose leaves correspond to vertices in $V union C$. For each edge $e$ in $T$, removing $e$ partitions the leaves into two sets $A$ and $B$. The _width_ of $T$ is:
  $
    "width"(T) = max_("edge" e "of" T) "cut-rank"(A_e, B_e)
  $
  
  The _rank width_ of $phi$ is:
  $
    "rw"(phi) = min_("rank decompositions" T) "width"(T)
  $
]

#theorem([Dynamic programming complexity])[
  The \#SAT problem for a formula with rank width $r$ can be solved in time $O(n dot 2^(3r))$ using dynamic programming along a rank decomposition of width $r$.
]

The rank width provides a finer measure of formula complexity compared to tree width. For example, formulas with dense clause structure can have bounded rank width but unbounded tree width.

==== Example: XOR constraints - Rank width vs Tree width

Consider a formula encoding XOR (parity) constraints. We express "$x_1 xor x_2 xor ... xor x_k = 0$" (even parity) in CNF, which requires exponentially many clauses but has nice algebraic structure.

For simplicity, consider a system of $n$ variables with $m$ XOR constraints, where each constraint involves $k$ variables. This can be represented as:
$
  phi = and.big_(i=1)^m "XOR-CNF"(S_i)
$
where $S_i subset {x_1, ..., x_n}$ with $|S_i| = k$.

*Example with $n=4$ variables, one global XOR constraint:*
$
  x_1 xor x_2 xor x_3 xor x_4 = 0 quad "(even parity)"
$

This is encoded in CNF as "an even number of variables must be true":
$
  phi = (overline(x)_1 or overline(x)_2 or overline(x)_3 or overline(x)_4) and 
        (x_1 or x_2 or overline(x)_3 or overline(x)_4) and ... quad "(8 clauses total)"
$

*Tree width analysis:* The incidence graph for XOR constraints is highly connected. Each of the 8 clauses connects to all 4 variables, creating a dense bipartite structure. The tree width is $Omega(k) = Omega(4) = 3$ or higher.

*Rank width analysis:* Consider the bipartition $A = {x_1, x_2}$ and $B = {x_3, x_4} union "clauses"$. 

The key insight: Over $bb(F)_2$, XOR constraints have *linear structure*. The effect of assigning $(x_1, x_2)$ on the clauses is determined by the *parity* $x_1 xor x_2$ alone!

The cut matrix $M_(A, B)$ has rows for ${x_1, x_2}$ and columns for ${x_3, x_4}$ and the clauses. Due to the XOR structure:
- All information about $(x_1, x_2)$ visible from $B$ is captured by their parity
- The rank is at most 2 (actually can be 1 for symmetric XOR formulas)

*General case:* For a system of $m$ linear constraints over $bb(F)_2$ on $n$ variables:
- Tree width: Can be $Theta(n)$ for dense constraint systems
- Rank width: At most $min(m, n)$, often much smaller due to linear dependencies

#figure(canvas(length: 1.5cm, {
  import draw: *
  
  // XOR constraint structure
  content((0, 2), text(11pt)[*XOR: $x_1 xor x_2 xor x_3 xor x_4 = 0$*])
  
  let var_y = 0.8
  let clause_y = 0
  
  // Variables
  for (i, x) in ((1, -2), (2, -0.7), (3, 0.7), (4, 2)) {
    circle((x, var_y), radius: 0.15, fill: blue.lighten(70%), name: "x"+str(i))
    content((x, var_y), text(8pt)[$x_#i$])
  }
  
  // Show 4 representative clauses (out of 8)
  for (i, x) in ((1, -1.8), (2, -0.6), (3, 0.6), (4, 1.8)) {
    circle((x, clause_y), radius: 0.12, fill: red.lighten(70%), name: "c"+str(i))
  }
  content((0, -0.5), text(7pt)[... 8 clauses ...])
  
  // Dense connections
  for i in range(1, 5) {
    for j in range(1, 5) {
      line("x"+str(i), "c"+str(j), stroke: (thickness: 0.3pt, paint: gray))
    }
  }
  
  // Cut
  line((0, -0.8), (0, 1.3), stroke: (thickness: 1.5pt, dash: "dashed", paint: red), name: "cut")
  content((0.3, 1.5), text(9pt, red)[cut])
  
  // Labels
  content((-1.3, 1.5), text(10pt, blue)[Side $A$])
  content((1.3, 1.5), text(10pt, green)[Side $B$])
  
  // Boundary info
  content((-1.3, -1.2), text(8pt)[4 assignments])
  content((-1.3, -1.5), text(8pt)[to $(x_1, x_2)$])
  content((0.5, -1.2), text(8pt, purple)[But only 2])
  content((0.5, -1.5), text(8pt, purple)[parities!])
  content((0.5, -1.8), text(7pt, purple)[$x_1 xor x_2 in {0,1}$])
}),
caption: [XOR constraint creates dense incidence graph (high tree width) but has low rank width because assignments are equivalent under parity over $bb(F)_2$.]
)

*Key insight:* For formulas encoding *linear algebraic structure* (systems of linear equations over $bb(F)_2$), rank width can be exponentially better than tree width:
- Tree width: $O(n)$ for dense systems
- Rank width: $O("rank of constraint matrix")$, often constant or logarithmic

Examples include error-correcting codes, cryptographic constraints, and parity games—all have natural $bb(F)_2$ structure that rank width captures perfectly.

==== Explicit computation: Rank width of a simple XOR constraint

Let's compute the rank width step-by-step for a minimal example: 3 variables with a single XOR constraint.

*Setup:* The constraint $x_1 xor x_2 xor x_3 = 0$ means "even number of variables are true".

$
  c_1 &: overline(x)_1 or overline(x)_2 or overline(x)_3 quad &"(0 true)"\
  c_2 &: x_1 or x_2 or overline(x)_3 quad &"(2 true: {x"_1",x"_2"})"\
  c_3 &: x_1 or overline(x)_2 or x_3 quad &"(2 true: {x"_1",x"_3"})"\
  c_4 &: overline(x)_1 or x_2 or x_3 quad &"(2 true: {x"_2",x"_3"})"
$

So $phi = c_1 and c_2 and c_3 and c_4$.

*Step 2: Incidence graph and matrix.* 

Variables: ${x_1, x_2, x_3}$, Clauses: ${c_1, c_2, c_3, c_4}$

The incidence matrix $M$ (variables as rows, clauses as columns):
$
  M = mat(
    1, 1, 1, 0;
    1, 1, 0, 1;
    1, 0, 1, 1
  ) quad "(over " bb(F)_2")"
$

For example, $M_(1,1) = 1$ because $x_1$ appears in $c_1$ (as $overline(x)_1$).

*Step 3: Bipartition $A = {x_1}$, $B = {x_2, x_3, c_1, c_2, c_3, c_4}$.*

The cut matrix $M_(A, B)$ (rows from $A$, columns from $B$):
$
  M_(A, B) = mat(augment: #2,
    1, 1, 1, 1, 1, 0
  )
$
(First 2 entries: edges to ${x_2, x_3}$; next 4: edges to clauses)

This is a single non-zero row:
$
  "cut-rank"({x_1}, B) = 1
$

*Step 4: Bipartition $A = {x_1, x_2}$, $B = {x_3, c_1, c_2, c_3, c_4}$.*

The cut matrix:
$
  M_(A, B) = mat(augment: #1,
    1, 1, 1, 1, 0;
    1, 1, 1, 0, 1
  )
$

Row reduce over $bb(F)_2$:
$
  mat(
    1, 1, 1, 1, 0;
    1, 1, 1, 0, 1
  ) arrow.r.long^("R2 + R1") mat(
    1, 1, 1, 1, 0;
    0, 0, 0, 1, 1
  )
$

Two linearly independent rows:
$
  "cut-rank"({x_1, x_2}, B) = 2
$

*Step 5: Check other bipartitions.*

By symmetry:
- Single variable partitions: cut-rank = 1
- Two variable partitions: cut-rank = 2  
- Partitions involving clauses: Similar analysis gives cut-rank ≤ 2

*Step 6: Conclusion.*

$
  "rw"(phi) = min_"all cuts" {"cut-rank"} = 1
$

The minimum is achieved by single-variable cuts.

*Why is this better than tree width?*

The incidence graph is a complete bipartite $K_(3,4)$ (each variable connects to most clauses). Any tree decomposition needs a large bag, giving tree width ≥ 2. But rank width = 1!

#figure(canvas(length: 1.5cm, {
  import draw: *
  
  // Left: tree decomposition
  content((-2, 2), text(11pt)[*Tree decomposition*])
  circle((-2, 1), radius: 0.5, name: "bag", fill: gray.lighten(80%))
  content((-2, 1), text(8pt)[$x_1,x_2,x_3$])
  content((-2, 0), text(9pt)[tw ≥ 2])
  
  // Right: rank decomposition  
  content((2, 2), text(11pt)[*Rank decomposition*])
  
  // Leaves
  circle((1, 0.5), radius: 0.15, fill: blue.lighten(70%), name: "x1")
  content((1, 0.5), text(8pt)[$x_1$])
  
  circle((2, 0.5), radius: 0.15, fill: blue.lighten(70%), name: "x23")
  content((2, 0.5), text(7pt)[$x_2$])
  
  circle((3, 0.5), radius: 0.15, fill: blue.lighten(70%), name: "x3")
  content((3, 0.5), text(7pt)[$x_3$])
  
  // Internal nodes
  circle((1.5, 1.3), radius: 0.1, fill: purple, name: "n1")
  circle((2.5, 1.3), radius: 0.1, fill: purple, name: "n2")
  circle((2, 2), radius: 0.1, fill: purple, name: "root")
  
  line("x1", "n1")
  line("x23", "n1")
  line("x3", "n2")
  line("n2", "n1", stroke: (dash: "dashed", thickness: 1.5pt, paint: red))
  
  content((1.7, 1.6), text(8pt, red)[rank=1])
  content((2, -0.2), text(9pt)[rw = 1])
}),
caption: [Tree width vs rank width for the XOR formula. Tree decomposition needs all variables in one bag (tw ≥ 2), but rank decomposition achieves width 1.]
)

*Key insight:* For the XOR constraint, two assignments $(x_1, x_2) = (0, 1)$ and $(1, 0)$ are equivalent because they have the same parity over $bb(F)_2$. The rank-1 structure captures this algebraic symmetry, collapsing 4 assignments into 2 equivalence classes.

#figure(canvas({
  import draw: *
  let d = 1.5
  
  // Draw a simpler incidence graph for XOR
  // Variables on left
  for (i, y) in ((1, 1.5), (2, 0.5), (3, -0.5), (4, -1.5)) {
    circle((0, y*0.6), radius: 0.2, name: "x"+str(i), fill: blue.lighten(70%))
    content((0, y*0.6), text(10pt)[$x_#i$])
  }
  
  // Representative clauses on right (showing structure)
  for (i, y) in ((1, 1.2), (2, 0.4), (3, -0.4), (4, -1.2)) {
    circle((2.5, y*0.6), radius: 0.15, name: "c"+str(i), fill: red.lighten(70%))
    content((2.5, y*0.6), text(8pt)[$c_#i$])
  }
  
  // Each clause connects to all 4 variables (dense)
  for i in range(1, 5) {
    for j in range(1, 5) {
      line("x"+str(i), "c"+str(j), stroke: (thickness: 0.4pt, paint: gray.lighten(20%)))
    }
  }
  
  content((1.25, -1.5), text(10pt)[Dense: Each clause])
  content((1.25, -1.8), text(10pt)[contains all 4 vars])
  
  content((1.25, 1.5), text(11pt)[Incidence graph])
}),
caption: [Incidence graph for XOR constraints. Each clause connects to many variables (dense structure), giving high tree width but low rank width due to $bb(F)_2$ linear structure.])

*Key insight:* For formulas encoding *linear algebraic structure* (systems of linear equations over $bb(F)_2$), rank width can be exponentially better than tree width:

=== Dynamic programming on rank decomposition for \#SAT

The key to exploiting small rank width is a dynamic programming algorithm that processes the rank decomposition tree bottom-up. The crucial observation is that we only need to track *equivalence classes* of partial assignments based on their *interaction with the other side of the cut*.

==== State representation

For each internal node of the rank decomposition corresponding to a bipartition $(A, B)$, we maintain a table of states. Each state is characterized by:
- A _boundary vector_ $bold(b) in bb(F)_2^r$ where $r = "cut-rank"(A, B)$
- The _count_ of satisfying partial assignments in subtree $A$ that produce boundary vector $bold(b)$

The boundary vector captures how the partial assignment in $A$ interacts with variables/clauses in $B$ through the linear span of the cut matrix $M_(A, B)$.

#definition([Boundary equivalence])[
  Two partial assignments $bold(x)_A, bold(x)'_A$ to variables in $A$ are _boundary equivalent_ if:
  $
    M_(A, B)^top bold(x)_A equiv M_(A, B)^top bold(x)'_A quad (mod 2)
  $
  where we treat the assignment as a vector in $bb(F)_2^(|A|)$.
]

*Intuition:* Imagine the bipartition $(A, B)$ as cutting the formula into two parts. To combine solutions from both sides later, we only need to know *how the assignment in $A$ affects clauses/variables in $B$*. 

The matrix-vector product $M_(A, B)^top bold(x)_A$ gives a vector in $bb(F)_2^(|B|)$ that encodes:
- For each element in $B$ (variable or clause), whether it's affected by the assignment in $A$

Two assignments $bold(x)_A$ and $bold(x)'_A$ are equivalent if they have the *same effect* on $B$. When combining with solutions from $B$ later, equivalent assignments behave identically—they satisfy the same clauses across the cut and interact with the same variables.

#figure(canvas(length: 1.8cm, {
  import draw: *
  
  // Draw the cut
  let h = 2
  line((0, -0.5*h), (0, 0.5*h), stroke: (thickness: 2pt, dash: "dashed", paint: red))
  
  // Left side (A)
  content((-1.5, 0.7*h), text(11pt, blue)[*Side A*])
  content((-1.5, 0.4*h), text(9pt)[$x_1 = 1, x_2 = 0$])
  content((-1.5, -0.1*h), text(9pt)[$x_1 = 0, x_2 = 1$])
  content((-1.5, -0.5*h), text(8pt, olive)[Different assignments])
  
  // Boundary vectors
  content((0, 0.4*h), text(9pt, red)[$arrow.r.long bold(b) = (1,0,1)$])
  content((0, -0.1*h), text(9pt, red)[$arrow.r.long bold(b) = (1,0,1)$])
  content((0, -0.8*h), text(8pt, purple)[Same boundary!])
  
  // Right side (B)
  content((1.5, 0.7*h), text(11pt, green)[*Side B*])
  content((1.5, 0.3*h), text(9pt)[sees: $(1,0,1)$])
  content((1.5, 0), text(8pt)[Can't distinguish])
  content((1.5, -0.3), text(8pt)[between $x_1=1,x_2=0$])
  content((1.5, -0.5), text(8pt)[and $x_1=0,x_2=1$])
}),
caption: [Boundary equivalence: Different assignments in $A$ that produce the same boundary vector are indistinguishable from the perspective of side $B$.]
)

*Example:* Consider variables $x_1, x_2 in A$ connected to clauses $c_1, c_2, c_3 in B$. Suppose:
- Clause $c_1$ contains both $x_1$ and $x_2$
- Clause $c_2$ contains only $x_1$
- Clause $c_3$ contains both $x_1$ and $x_2$

In $bb(F)_2$ (mod 2 arithmetic):
- Assignment $(x_1=1, x_2=0)$ affects $c_1$ by 1, $c_2$ by 1, $c_3$ by 1 → boundary $(1,1,1)$
- Assignment $(x_1=0, x_2=1)$ affects $c_1$ by 1, $c_2$ by 0, $c_3$ by 1 → boundary $(1,0,1)$
- Assignment $(x_1=1, x_2=1)$ affects $c_1$ by $1+1=0$ (mod 2), $c_2$ by 1, $c_3$ by 0 → boundary $(0,1,0)$

Assignments with the same boundary are *interchangeable* when combining with solutions from $B$—they contribute identically to clause satisfaction across the cut.

Since the cut-rank is $r$, there are at most $2^r$ distinct boundary vectors, regardless of the size of $A$. This is the key to efficiency!

==== Basis transformation perspective

Boundary equivalence can be elegantly interpreted as a *change of basis* for the variables in $A$.

*Original representation:* Variables $bold(x)_A in bb(F)_2^(|A|)$ (standard basis)
- There are $2^(|A|)$ possible assignments
- Each assignment is a vertex of the $|A|$-dimensional Boolean hypercube

*Transformed representation:* Boundary vectors $bold(b) in bb(F)_2^(|B|)$ via $bold(b) = M_(A,B)^top bold(x)_A$
- The image lies in a rank-$r$ subspace of $bb(F)_2^(|B|)$
- There are at most $2^r$ distinct boundary vectors

The linear map $M_(A,B)^top : bb(F)_2^(|A|) arrow.r bb(F)_2^(|B|)$ projects the $|A|$-dimensional space onto its $r$-dimensional image. We can decompose:
$
  bb(F)_2^(|A|) = "ker"(M_(A,B)^top) plus.o.big "complementary space"
$
where $dim("ker"(M_(A,B)^top)) = |A| - r$.

*What this means:*
- Variables in $A$ can be split into two groups:
  1. *$r$ "boundary-relevant" degrees of freedom* that affect side $B$
  2. *$|A| - r$ "internal" degrees of freedom* invisible to side $B$
  
- We can choose a new basis ${bold(y)_1, ..., bold(y)_r, bold(z)_1, ..., bold(z)_(|A|-r)}$ where:
  - ${bold(y)_1, ..., bold(y)_r}$ span the image of $M_(A,B)^top$ (the "boundary basis")
  - ${bold(z)_1, ..., bold(z)_(|A|-r)}$ span the kernel (the "internal basis")

*Dynamic programming in the new basis:*
- State space: Only track the $r$ boundary variables ${bold(y)_1, ..., bold(y)_r}$ → $2^r$ states
- Internal variables: Sum over all $2^(|A|-r)$ configurations of ${bold(z)_1, ..., bold(z)_(|A|-r)}$ for each boundary state

#figure(canvas(length: 1.5cm, {
  import draw: *
  
  // Left side: original space
  content((-2, 2.5), text(11pt)[*Original variables*])
  content((-2, 2), text(10pt)[$bold(x)_A in bb(F)_2^(|A|)$])
  
  // Draw hypercube representation
  let s = 0.4
  for i in range(0, 4) {
    for j in range(0, 4) {
      circle((-2.5 + i*s, 0.5 + j*s), radius: 0.05, fill: blue.lighten(70%))
    }
  }
  content((-2, -0.3), text(9pt)[$2^(|A|)$ assignments])
  
  // Arrow
  content((0, 1), text(12pt)[$M_(A,B)^top$])
  line((-1, 1), (1, 1), mark: (end: ">"))
  
  // Right side: transformed space  
  content((2, 2.5), text(11pt)[*Boundary vectors*])
  content((2, 2), text(10pt)[$bold(b) in bb(F)_2^r$])
  
  // Draw reduced space
  for i in range(0, 2) {
    for j in range(0, 2) {
      circle((1.5 + i*0.6, 0.5 + j*0.6), radius: 0.08, fill: red.lighten(60%), stroke: red)
    }
  }
  content((2, -0.3), text(9pt)[$2^r$ equivalence classes])
  
  // Labels
  content((-2, -1), text(9pt, blue)[High dimension])
  content((2, -1), text(9pt, red)[Low dimension])
}),
caption: [Boundary equivalence as dimensionality reduction: The linear map $M_(A,B)^top$ projects $2^(|A|)$ assignments onto $2^r$ equivalence classes, where $r = "rank"(M_(A,B))$.]
)

*XOR example revisited:* For $x_1 xor x_2 xor x_3 = 0$ with partition $A = {x_1, x_2}$:
- Original basis: $(x_1, x_2) in bb(F)_2^2$ → 4 assignments
- Boundary basis: $(x_1 xor x_2) in bb(F)_2^1$ → 2 equivalence classes
- The parity $y_1 = x_1 xor x_2$ is the single "boundary-relevant" variable
- We can choose $z_1 = x_1$ as the "internal" variable (invisible from side $B$)
- DP state space: Only track $y_1$ (2 states), sum over $z_1$ internally

This basis transformation perspective explains why rank width is so powerful: it automatically discovers the *minimal set of effective degrees of freedom* at the boundary, exploiting algebraic structure that tree-based methods cannot see.

==== Connection to tensor networks: SVD compression

The basis transformation view of boundary equivalence has a direct correspondence with *SVD compression* in tensor networks (MPS/TTN)!

*Encoding \#SAT as a tensor network:*

Each variable $x_i$ and clause $c_j$ can be encoded as a tensor:
- Variable tensor: $T^(x_i)_(s_i) = 1$ for $s_i in {0, 1}$ (identity, dimension 2)
- Clause tensor: $C^(c_j)_(s_(i_1), ..., s_(i_k)) = cases(1 "if clause satisfied", 0 "otherwise")$

The model count is the tensor network contraction:
$
  \#"SAT"(phi) = sum_(s_1, ..., s_n in {0,1}) product_j C^(c_j)(...) = "contract"(T N)
$

*Rank decomposition ↔ Tensor tree network:*

A rank decomposition naturally defines a binary tree tensor network (TTN):
- Leaves: Variable/clause tensors
- Internal nodes: Partial contractions over subtrees
- Each edge: A bond connecting two subtrees

*Boundary equivalence ↔ SVD at bonds:*

At each cut $(A, B)$ in the rank decomposition:

#table(
  columns: (auto, auto),
  align: left,
  [*Rank decomposition*], [*Tensor network*],
  [Assignment $bold(x)_A in {0,1}^(|A|)$], [Physical indices in subtree $A$],
  [Boundary vector $bold(b) in bb(F)_2^r$], [Virtual bond index (dimension $2^r$)],
  [Linear map $M_(A,B)^top$], [Unfolding + SVD of contracted tensor],
  [Cut-rank $r$], [Bond dimension after compression],
  [Basis transformation], [Gauge transformation to canonical form]
)

*The SVD connection:*

Let $T_A$ be the tensor representing all configurations in subtree $A$, with shape $2^(|A|) times dots.c$ (physical indices on one side, rest on the other). When we contract towards the cut:

1. *Unfold* $T_A$ into a matrix with rows = internal indices, columns = boundary indices
2. *Compute SVD*: $T_A = U Sigma V^dagger$
3. *Truncate* to rank $r$: Keep only singular values corresponding to non-zero entries
4. The rank $r$ equals the cut-rank over $bb(F)_2$ for Boolean tensors!

The boundary vector $bold(b)$ corresponds to the *virtual bond index* after SVD compression.

#figure(canvas(length: 1.4cm, {
  import draw: *
  
  // Left subtree (side A)
  content((-2.5, 2), text(10pt)[*Subtree A*])
  for (i, y) in ((1, 0.5), (2, 0), (3, -0.5)) {
    circle((-3, y), radius: 0.15, fill: blue.lighten(70%), name: "a"+str(i))
    content((-3, y), text(7pt)[$x_#i$])
  }
  
  // Contracted tensor
  circle((-1.5, 0), radius: 0.3, fill: purple.lighten(70%), name: "TA")
  content((-1.5, 0), text(9pt)[$T_A$])
  
  for i in range(1, 4) {
    line("a"+str(i), "TA", stroke: (thickness: 0.5pt))
  }
  
  // Bond (before compression)
  line((-1.2, 0), (-0.3, 0), stroke: (thickness: 3pt, paint: gray), name: "bond1")
  content((-0.75, 0.4), text(8pt)[dim $2^(|A|)$])
  
  // SVD arrow
  content((0.5, 0), text(11pt)[SVD])
  line((1, 0.2), (1.5, 0.2), mark: (end: ">"))
  
  // After compression
  circle((2.2, 0), radius: 0.3, fill: purple.lighten(70%), name: "U")
  content((2.2, 0), text(9pt)[$U$])
  line((2.5, 0), (3.2, 0), stroke: (thickness: 1.5pt, paint: red), name: "bond2")
  content((2.85, 0.4), text(8pt, red)[dim $2^r$])
  
  content((3.5, 0.4), text(9pt)[$dots.c arrow.r$])
  content((2.2, -0.8), text(9pt)[Compressed])
  
  content((-1.5, -1.5), text(9pt)[Basis: original])
  content((2.2, -1.5), text(9pt)[Basis: boundary])
}),
caption: [Tensor network perspective: Contracting subtree $A$ produces tensor $T_A$. SVD at the cut compresses bond dimension from $2^(|A|)$ to $2^r$, corresponding to the basis transformation from original variables to boundary equivalence classes.]
)

*Why Boolean tensors have special structure:*

For \#SAT, tensors have entries in ${0, 1}$. Over $bb(F)_2$:
- Many configurations may map to the same boundary state due to linear dependencies
- The effective rank is often much smaller than the naive dimension $2^(|A|)$
- SVD automatically finds this compressed representation

*MPS example:* For the XOR constraint $x_1 xor x_2 xor x_3 = 0$, with partition between $(x_1, x_2)$ and $(x_3)$:
- Variables in left part: ${x_1, x_2}$
- Naive bond dimension (TN): $2^2 = 4$ configurations of $(x_1, x_2)$
- But rank width considers the incidence graph including clauses!
- Rank width at this cut = 1 (computed from the full incidence matrix)
- DP state space: $2^1 = 2$ states (corresponding to parity equivalence)
- The DP state space size relates to, but is not exactly, the tensor bond dimension

*Important distinction:* Rank width and bond dimension are related but NOT the same!

#proposition([Rank width bounds DP state space])[
  For the rank-decomposition-based DP algorithm for \#SAT:
  - The rank width $r$ at a cut $(A, B)$ bounds the number of _equivalence classes_ of partial assignments
  - The DP state space has dimension $O(2^r)$ 
  - This is NOT the same as tensor network bond dimension, which lives in configuration space
]

*The key distinction:*

#table(
  columns: (auto, auto),
  align: left,
  [*Rank width $r$*], [*Bond dimension $chi$*],
  [Computed from incidence matrix $M_(A,B)$], [Dimension of virtual index in tensor],
  [Rank over $bb(F)_2$ of graph structure], [Size of index in configuration space],
  [Measures linear dependencies in constraints], [Measures entanglement/correlation],
  [Partition includes variables AND clauses], [Partition of physical variables only],
  [DP state space: $O(2^r)$], [Tensor bond: $O(2^(|A_"var"|))$ before compression],
)

*What's actually happening:*

When we partition the incidence graph into $(A, B)$:
- $A$ contains some variables ${x_1, ..., x_k}$ and some clauses ${c_1, ..., c_m}$
- Rank width $r = "rank"(M_(A,B))$ over $bb(F)_2$
- DP tracks equivalence classes based on how assignments affect the other side
- Number of DP states $approx 2^r$, but this is about the _information flow_, not directly bond dimension

In tensor network language:
- Bond dimension = size of virtual index between contracted parts
- For variables in $A$: naive bond dimension is $2^k$ (all configurations)
- Compression reduces this, but the achievable compression depends on:
  1. How many variables are in the cut
  2. The constraint structure (captured partially by rank width)
  
The relationship is subtle: low rank width implies the DP is efficient, which means the _effective_ information crossing the boundary is small, but this doesn't directly translate to a simple formula for tensor bond dimension.

This connection reveals that:
- Rank decomposition DP = Tensor tree network contraction with compression
- Low rank width = DP state space is small ($O(2^r)$ equivalence classes)
- The $bb(F)_2$ structure in Boolean problems enables efficient state space compression

==== Why rank width matters: Finding the optimal tree structure

*Question:* If SVD can capture low-rankness at any bond, why do we need rank width? Can't we just apply SVD to an arbitrary tensor tree?

*Answer:* Rank width solves a fundamentally different problem: _finding the optimal tree structure itself_.

The key distinction:

#table(
  columns: (auto, auto),
  align: left,
  [*SVD (Computational tool)*], [*Rank width (Structural parameter)*],
  [Compresses a *given* bond], [Finds the *best tree decomposition*],
  [Local operation at one cut], [Global optimization over all trees],
  [Answers: "How low-rank is *this* bond?"], [Answers: "What tree gives the *lowest* max rank?"],
  [Algorithmic technique], [Complexity measure],
)

*The problem with arbitrary trees:*

Consider the XOR constraint $x_1 xor x_2 xor x_3 = 0$:

*Bad tree structure (sequential):*
```
    x₁ ─ x₂ ─ x₃ ─ clauses
```
- Cut after $x_1$: rank = 1 ✓
- Cut after $x_2$: Must track $(x_1, x_2)$ interactions with $(x_3, "clauses")$ → rank = 2
- Maximum rank = 2

*Good tree structure (balanced):*
```
      ┌─x₁
    ┌─┤
  ──┤ └─x₂
    └─x₃ + clauses
```
- Every cut: rank = 1 ✓
- Maximum rank = 1

SVD at each bond gives the same local compressions, but the *global* complexity depends on the *tree structure*!

*Rank width = Optimal tree structure:*

Rank width is the solution to the optimization problem:
$
  "rw"(phi) = min_("binary trees" T) max_("cuts" e "in" T) "rank"(M_(A_e, B_e))
$

This is a _global optimization_ over exponentially many possible tree structures. Finding the optimal tree is itself a hard problem (but tractable for bounded rank width using dynamic programming on the powerset).

*Why this matters for complexity:*

#table(
  columns: (auto, auto, auto),
  align: left,
  [*Tree structure*], [*Max bond dim*], [*Complexity*],
  [Arbitrary tree], [$2^(O(n))$ worst case], [Exponential],
  [Tree width decomposition], [$2^("tw")$], [Exponential in tw],
  [Optimal rank decomp], [$2^("rw")$], [Exponential in rw],
)

For XOR formulas:
- Arbitrary tree: max bond dim = $O(n)$ → complexity $2^(O(n))$
- Optimal rank decomp: max bond dim = $O(1)$ → complexity $O(n)$

*Practical implications:*

1. *Without rank width analysis:* 
   - Use heuristic tree structure
   - Apply SVD at each bond
   - No guarantees on overall complexity
   - May still be exponential

2. *With rank width analysis:*
   - Find provably optimal tree structure
   - Guarantees $max "bond" = 2^("rw")$
   - Polynomial time for constant rank width
   - Can prove tractability a priori

*The deeper insight:*

Rank width is not just about compression—it's about _problem structure_. It tells us:
- Whether a problem has hidden algebraic structure
- What tree decomposition exploits that structure optimally  
- What complexity class the problem belongs to

SVD is the *mechanism* for exploiting low-rankness locally. Rank width is the *theory* that tells us:
1. Whether low-rankness exists globally
2. How to find the tree structure that achieves it
3. What complexity to expect

#theorem([Rank width as complexity measure])[
  A problem with rank width $r$ can be solved in time $O(n times 2^(O(r^2)))$. This is a _structural_ complexity bound, independent of the specific algorithm used—any method that exploits the optimal tree structure achieves this.
]

In summary: _SVD is the tool, rank width is the blueprint that tells you where and how to use it optimally._

==== Algorithm outline

*Step 1: Initialize leaves.* 
- For variable leaf $x_i$: Create states for $x_i = 0$ and $x_i = 1$, each with count = 1
- For clause leaf $c_j$: Create one state (clauses don't get assigned), count = 1

*Step 2: Merge internal nodes.* At each internal node combining subtrees with bipartition $(A_"left", A_"right")$ into $A$:
- For each pair of boundary vectors $(bold(b)_"left", bold(b)_"right")$:
  + Compute the combined boundary vector $bold(b)_"new"$ for the parent cut
  + Check constraints between left and right (clauses satisfied, etc.)
  + If valid: Update $C[bold(b)_"new"] arrow.l C[bold(b)_"new"] + C_"left"[bold(b)_"left"] dot C_"right"[bold(b)_"right"]$

*Step 3: Root node.* At the root, the boundary is empty (no cut). Sum counts over all boundary vectors that represent globally satisfying assignments.

==== Complexity analysis

For a rank decomposition of width $r$:
- Number of states per node: $O(2^r)$
- Number of nodes: $O(n + m)$ where $n$ is variables and $m$ is clauses
- Merge operation: $O(2^(2r))$ to combine two subtrees
- Total time: $O((n + m) dot 2^(2r))$

For formulas with constant rank width, this gives polynomial time!

#theorem([Efficiency of rank-decomposition DP])[
  For \#SAT on a formula with rank width $r$:
  - States per node: $2^r$ (independent of subtree size)
  - Total complexity: $O(n dot 2^(3r))$
  - For constant rank width, this is polynomial time
]

The exponential improvement comes from the fact that rank width captures the _information bottleneck_ between subtrees more precisely than tree width. The rank over $bb(F)_2$ identifies which partial assignments are equivalent in terms of their effect on satisfying clauses across the cut.
