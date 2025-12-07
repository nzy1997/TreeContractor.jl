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

The rank width provides a finer measure of formula complexity compared to tree width. For example, formulas with dense clause structure can have bounded rank width but unbounded tree width.

==== Example: Parity constraint (XOR)

Consider the constraint $x_1 xor x_2 xor x_3 = 0$ (even parity). In CNF, this requires 4 clauses, each connecting all 3 variables (e.g., $x_1 or x_2 or overline(x)_3$).

*Tree Width (Structural Complexity):*
The incidence graph is a complete bipartite graph $K_(3,4)$ because every variable connects to every clause.
- *Structure:* A dense web of connections.
- *Implication:* Any tree decomposition must put all 3 variables in a single bag.
- *Tree Width:* $>= 2$.

*Rank Width (Algebraic Complexity):*
We analyze the information flow across a cut over $bb(F)_2$. Consider the partition $A={x_1}$ and $B={x_2, x_3, "clauses"}$.
- *Effect of $A$:* The value of $x_1$ only sends 1 bit of information to the rest of the system: "what is my contribution to the parity?".
- *Cut Matrix:* In the incidence matrix, the row for $x_1$ represents its connections to clauses. Despite connecting to all 4 clauses, this row vector lives in a 1-dimensional subspace relative to the rest of the graph's structure (due to linear dependencies over $bb(F)_2$).
- *Rank Width:* 1.

*Key Takeaway:* Tree width sees "many edges" (high complexity). Rank width sees "linear dependency" (low complexity). This makes rank width powerful for problems with hidden algebraic structure like XOR-SAT.

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
caption: [Tree width vs rank width for the XOR formula. Tree decomposition needs all variables in one bag (tw ≥ 2) due to dense connections, but rank decomposition achieves width 1 by exploiting linear dependency.]
)

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

==== Connection to tensor networks and SVD

We can now rigorously connect three concepts: rank decomposition, tensor networks, and SVD compression.

*1. Encoding \#SAT as a tensor network*

Each variable $x_i$ and clause $c_j$ becomes a tensor in a network:
- Variable tensor $T^(x_i)$: Rank-1 tensor (vector) of dimension 2, representing physical state.
- Clause tensor $C^(c_j)$: Rank-$k$ tensor where $k$ is clause width. $C = 1$ if satisfied, 0 otherwise.

The model count is the contraction of this entire network:
$
  \#"SAT"(phi) = "contract"(product_i T^(x_i) times product_j C^(c_j))
$

*2. Tree structure induces contraction order*

A rank decomposition tree $T$ defines a *hierarchical contraction order*:
- Leaves are individual tensors.
- Internal nodes represent partial contractions of subtrees.
- Edges represent virtual bonds between intermediate tensors.

*3. Rank width bounds the bond dimension*

At any cut $(A, B)$ in the tree, we have a partially contracted tensor $cal(T)_A$ representing the left subtree.

#proposition([Equivalence of Rank Width and Compressed Bond Dimension])[
  Let $r = "rank"(M_(A,B))$ be the cut-rank of the incidence matrix over $bb(F)_2$.
  Then, the tensor $cal(T)_A$ can be compressed via SVD to a bond dimension $chi = 2^r$.
]

*Proof Sketch:*
The incidence matrix $M_(A,B)$ captures the *linear dependencies* in how variables in $A$ interact with clauses in $B$.
- If $M_(A,B)$ has rank $r$, there are $2^r$ equivalence classes of assignments in $A$.
- These equivalence classes form an orthogonal basis for the interaction space.
- SVD on the tensor $cal(T)_A$ (unfolded) will find exactly $2^r$ non-zero singular values corresponding to these classes.

*Example: XOR constraint $x_1 xor x_2 xor x_3 = 0$*

Let's compare the naive bond dimension vs. the compressed one for the cut $A={x_1, x_2}, B={x_3, "clauses"}$.

*A. Naive Contraction (No SVD)*
- Tensor $cal(T)_A$ carries indices for $x_1, x_2$.
- Bond dimension $chi_"naive" = 2^2 = 4$.
- States: $(0,0), (0,1), (1,0), (1,1)$.

*B. SVD Compression (Tensor Network)*
- Unfold $cal(T)_A$ to matrix.
- SVD reveals only 2 non-zero singular values.
- Compressed bond dimension $chi_"SVD" = 2$.
- The relevant states are parity classes: $x_1 xor x_2 = 0$ and $x_1 xor x_2 = 1$.

*C. Rank Width Analysis (Decision Diagram)*
- For this cut $A={x_1, x_2}$, the incidence matrix rank is $r = 2$.
- DP state space bound: $2^r = 4$.
- *Note:* Rank width purely based on graph structure is an upper bound. SVD can be tighter (2) by exploiting value symmetries.
- *However*, rank width guides us to find a *better tree* (e.g., cutting at $x_1$) where $r=1$, matching the optimal SVD complexity!

*4. Why Rank Width? (Global Optimization)*

If SVD can compress bonds, why do we need rank width?

*Answer:* *Global Structure Optimization.*

SVD is a *local* tool—it compresses a specific bond given a specific tree. Rank width is the *global* measure that tells us *which tree* allows for the best compression everywhere.

#table(
  columns: (auto, auto),
  align: left,
  [*Local Tool (SVD)*], [*Global Theory (Rank Width)*],
  [Given a tree, optimizes bond $e$], [Finds the optimal tree $T$],
  [Computational reduction], [Structural complexity bound],
  [$chi_e arrow.r 2^(r_e)$], [$min_T max_e r_e$],
)

For the XOR example:
- *Bad Tree:* Linear chain $x_1 - x_2 - x_3$. Max cut-rank = 2. Max bond dim = 4.
- *Good Tree:* Balanced decomposition. Max cut-rank = 1. Max bond dim = 2.

Rank width theory guarantees that a "good tree" exists and bounds its worst-case complexity. Tensor networks with SVD provide the practical machinery to execute the contraction efficiently on that tree.

#theorem([Unified View])[
  The rank-decomposition DP algorithm is mathematically equivalent to contracting a Tensor Tree Network (TTN) where every bond is compressed to its optimal rank via SVD. The "rank width" of the formula corresponds to the logarithm of the maximum bond dimension required in the optimal TTN.
]

==== Algorithm outline

The dynamic programming proceeds bottom-up on the rank decomposition tree:

*Step 1: Initialize leaves.* 
- For variable leaf $x_i$: Create states for $x_i = 0$ and $x_i = 1$, each with count = 1.
- For clause leaf $c_j$: Create one state (clauses don't get assigned), count = 1.

*Step 2: Merge internal nodes.* 
At a node combining subtrees $L$ and $R$:
- Iterate over all pairs of states $(b_L, b_R)$ from children.
- Compute new boundary $b_"new" = b_L + b_R$ (or appropriate linear combination).
- Check consistency (e.g., if any fully-formed clauses are satisfied).
- Update: $C(b_"new") arrow.l C(b_"new") + C_L(b_L) dot C_R(b_R)$.

*Step 3: Root node.* 
Sum counts over all valid states at the root (usually the empty boundary vector $bold(0)$).

==== Complexity analysis

For a rank decomposition of width $r$:
- *State space*: Each node stores at most $2^r$ states.
- *Merge cost*: Combining two tables requires iterating pairs or matrix multiplication, taking $O(2^(3r))$ time (or $O(2^(2r))$ with sparse interactions).
- *Total time*: With $N$ nodes in the decomposition, the complexity is $O(N dot 2^(3r))$.

This matches the complexity of contracting the corresponding tensor network with bond dimension $2^r$.

#theorem([Efficiency of Rank Width DP])[
  The \#SAT problem for a formula with rank width $r$ can be solved in time $O(n dot 2^(3r))$. This is polynomial for constant $r$, providing a tractable approach for structured instances like XOR formulas where tree width fails.
]
