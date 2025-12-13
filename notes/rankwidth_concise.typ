#import "@preview/cetz:0.4.0": canvas, draw, tree, coordinate
#import "@preview/cetz-plot:0.1.2": *
#import "@preview/ctheorems:1.1.3": *
#import "@preview/algorithmic:1.0.7": *

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


== \#SAT and rank width

The \#SAT (model counting) problem asks: given a Boolean formula $phi$ in CNF, how many satisfying assignments exist? Model counting is \#P-complete, but becomes tractable for formulas with bounded rank width.

=== Rank width

#definition([Incidence graph])[
  For a CNF formula $phi$ with variables $V$ and clauses $C$, the _incidence graph_ $G = (V union C, E)$ is a bipartite graph where each variable $v in V$ is connected to each clause $c in C$ that contains $v$ or $overline(v)$.
  
  *Note:* The incidence graph captures *syntactic connectivity* (which variables appear in which clauses) rather than *semantic polarity* (whether literals are positive or negative). This is intentional because:
  1. Rank width measures structural complexity of variable-clause interactions
  2. Over $bb(F)_2$, both $x$ and $overline(x)$ represent the same variable contributing to constraints
  3. Polarity affects constraint satisfaction checking (handled separately in DP), not the cut-rank computation
]

#definition([Rank width])[
  Let $G = (V union C, E)$ be the incidence graph. For a bipartition $(A, B)$ of vertices, the _cut-rank_ is $"cut-rank"(A, B) = "rank"_(bb(F)_2)(M_(A, B))$, where $M_(A, B)$ is the submatrix with rows indexed by $A$ and columns by $B$, computed over $bb(F)_2$.
  
  A _rank decomposition_ is a binary tree $T$ with leaves in $V union C$. The _width_ is $max_("edge" e) "cut-rank"(A_e, B_e)$. The _rank width_ $"rw"(phi)$ is the minimum width over all rank decompositions.
]

*Why $bb(F)_2$?* Operations over the binary field directly match Boolean logic: $x + y equiv x xor y$ and $x dot y equiv x and y$. Rank over $bb(F)_2$ captures linear dependencies in Boolean constraints (e.g., XOR is addition mod 2), revealing algebraic structure invisible to tree-based methods.

*Example - Polarity doesn't affect incidence structure:*
Consider two clauses:
- $c_1 = x_1 or x_2 or overline(x)_3$ (mixed polarities)
- $c_2 = overline(x)_1 or overline(x)_2 or x_3$ (different polarities)

Both clauses connect to the same three variables ${x_1, x_2, x_3}$ in the incidence graph. The cut matrix entry $M_(i,j) = 1$ simply indicates "$x_i$ appears in $c_j$" (regardless of sign).

*Where polarity matters:*
During DP, when checking if a partial assignment satisfies a clause, we use the actual polarities:
- For $c_1$: need $(x_1=1) or (x_2=1) or (x_3=0)$
- For $c_2$: need $(x_1=0) or (x_2=0) or (x_3=1)$

But for rank-width computation (structure analysis), we only care that both clauses "involve" variables ${x_1, x_2, x_3}$.

=== Example: System of XOR constraints

*Problem:* Count satisfying assignments for a system of XOR equations over 5 variables:
$
  cases(
    x_1 xor x_2 xor x_3 = 0,
    x_2 xor x_3 xor x_4 = 0,
    x_3 xor x_4 xor x_5 = 0
  )
$

Each equation becomes 4 CNF clauses, giving 12 clauses total. The incidence graph is dense with 60 edges.

*Known answer:* 8 satisfying assignments (solutions to the linear system over $bb(F)_2$).

*Tree width analysis:* The incidence graph has cliques and dense connections. Any tree decomposition requires width $>= 3$, giving $O(2^3) = O(8)$ states at separators.

*Rank width analysis:* Consider partition $A = {x_1, x_2}$, $B = {x_3, x_4, x_5, "clauses"}$. 

The cut matrix captures which clauses in $B$ connect to variables in $A$. Each equation becomes 4 CNF clauses, so we have 12 clauses: $c_1, dots, c_12$.

*Variable appearances:*
- $x_1$: appears in $c_1, c_2, c_3, c_4$ (equation 1 only)
- $x_2$: appears in $c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8$ (equations 1 and 2)

The complete cut matrix is:
$
  M_(A, B) = mat(
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0;
  ) in bb(F)_2^(2 times 12)
$
where columns correspond to clauses $c_1, c_2, dots, c_12$ respectively.

Over $bb(F)_2$, the second row is NOT a multiple of the first:
$
  "row"_2 = "row"_1 + (0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0)
$
Thus the rows are linearly independent, giving:
$
  "rank"_(bb(F)_2)(M_(A, B)) = 2
$

*Key insight:* Despite $x_1, x_2$ appearing in 8 of 12 clauses (dense connectivity), they only send 2 bits of information: their individual parities. The rank width $<= 2$, giving $O(2^2) = O(4)$ states - a 2× reduction from tree width.

*Scaling:* For a chain of $n$ variables with $n-2$ overlapping XOR constraints:
- Naive: $O(2^n)$ 
- Tree width: $O(n dot 2^3)$ (width 3 from overlapping constraints)
- Rank width: $O(n dot 2^2)$ (width 2 from parity structure)

For $n=10$: Tree-based DP uses $~10,000$ states, rank-based uses $~40$ states (250× speedup).

==== Computation trace (sketch)

*Phase 1* (leaves): Each variable: 2 states $(0,1), (1,1)$.

*Phase 2* (merge $x_1, x_2$): Boundary vector $bold(b) = (x_1, x_2) in bb(F)_2^2$:
- 4 boundary states: $(0,0), (0,1), (1,0), (1,1)$ each with count 1

*Phase 3* (merge with $x_3$): Combined with constraint $x_1 xor x_2 xor x_3 = 0$:
- Filter valid combinations based on first constraint
- New boundary: $(x_2, x_3)$ for next constraint $x_2 xor x_3 xor x_4 = 0$
- State space remains $2^2 = 4$

*Phase 4-5* (Continue merging): Each merge maintains 4 boundary states due to rank 2 structure.

*Result:* Final count = 8 ✓

==== Complexity comparison

#table(
  columns: (auto, auto, auto),
  align: left,
  [*Approach*], [*State Space (5 vars)*], [*Scaling for $n$ variables*],
  [Naive], [$2^5 = 32$], [$O(2^n)$],
  [Tree-based (tw=3)], [$O(8)$], [$O(n dot 2^3)$],
  [Rank-based (rw=2)], [$O(4)$], [$O(n dot 2^2)$],
)

*Concrete advantage:* For $n=10$ variables with overlapping XOR constraints:
- Tree-based: $~10 dot 2^3 = 80$ node-states
- Rank-based: $~10 dot 2^2 = 40$ node-states
- For $n=20$: Tree-based $~160$ vs Rank-based $~80$ (advantage grows linearly)

*Key:* Tree width counts structural connectivity. Rank width counts algebraic degrees of freedom over $bb(F)_2$. For systems of linear equations over $bb(F)_2$, rank width captures the true constraint complexity, independent of how densely the CNF encoding connects variables.

=== Algorithm: Rank-based model counting

#import "@preview/algorithmic:1.0.7"
#algorithm({
  import algorithmic: *
  Function("RankDP", ([$v$], [$phi$]), {
    Comment([Recursively compute boundary states for subtree rooted at $v$])
    
    If([$v$ is variable leaf $x_i$], {
      Return[${(0, 1), (1, 1)}$ #CommentInline([Two values, count 1 each])]
    })
    
    If([$v$ is clause leaf $c_j$], {
      Return[${(emptyset, 1)}$ #CommentInline([Clause carries no initial state])]
    })
    
    Comment([Internal node: recursively process children])
    Assign([$"States"_L$], [RankDP(left($v$), $phi$)])
    Assign([$"States"_R$], [RankDP(right($v$), $phi$)])
    
    Comment([Compute cut-rank for partition at $v$])
    Assign([$(A, B)$], [partition induced by removing $v$ from tree])
    Comment([Build cut matrix: rows = vertices in $A$, cols = vertices in $B$])
    Assign([$M_(A,B)[i,j]$], [1 if edge $(i,j)$ crosses cut, 0 otherwise])
    Assign([$r$], [$"rank"_(bb(F)_2)(M_(A,B))$ #CommentInline([At most $2^r$ states])])
    
    Assign([$"States"_v$], [$emptyset$])
    For([$(bold(b)_L, c_L) in "States"_L$, $(bold(b)_R, c_R) in "States"_R$], {
      Comment([Merge boundary vectors])
      Assign([$bold(b)$], [$phi.alt(bold(b)_L, bold(b)_R)$ #CommentInline([Problem-specific])])
      
      If([constraints on $A$ satisfied], {
        Assign([$"States"_v[bold(b)]$], [$"States"_v[bold(b)] + c_L times c_R$])
      })
    })
    
    Return[$"States"_v$]
  })
  
  Assign([$"root_states"$], [RankDP(root($T$), $phi$)])
  Return[$sum_("valid" bold(b)) "root_states"[bold(b)]$]
})

*Key operations:*
- *Base cases*: Variable leaves return 2 states (value + count), clause leaves return 1 state
- *Recursive merge*: Process left and right subtrees first (post-order traversal)
- *Cut matrix $M_(A,B)$*: At each internal node $v$, construct the adjacency matrix between partition sets $A$ and $B$. The rank $r = "rank"_(bb(F)_2)(M_(A,B))$ bounds the dimension of boundary state space
- *Boundary vectors $bold(b) in bb(F)_2^r$*: For partition $(A, B)$ with cut matrix $M_(A,B)$, the boundary vector is $bold(b) = M_(A,B)^top bold(x)_A$ where $bold(x)_A in {0,1}^(|A|)$ is the assignment to variables in $A$, and $r = "rank"_(bb(F)_2)(M_(A,B))$. This projects the $2^(|A|)$ possible assignments onto $2^r$ equivalence classes
- *Boundary merge $phi.alt(bold(b)_L, bold(b)_R)$*: Combines boundary vectors from left and right subtrees:
  - Input: $bold(b)_L in bb(F)_2^(r_L)$ where $r_L$ is the rank of left subtree's cut
  - Input: $bold(b)_R in bb(F)_2^(r_R)$ where $r_R$ is the rank of right subtree's cut
  - Expand: Convert $bold(b)_L, bold(b)_R$ back to partial assignments $bold(x)_L, bold(x)_R$
  - Check: All constraints *internal* to merged subtree satisfied (clauses connecting $L$ and $R$)
  - Project: Compute new boundary $bold(b) = M_(A^',B^')^top (bold(x)_L union bold(x)_R)$ for parent cut $(A^', B^')$
  - Output: $bold(b) in bb(F)_2^r$ where $r = "rank"(M_(A^',B^'))$
  - Example (XOR): For partition $A = {x_1, x_2}$, $B = {x_3, "clauses"}$ with $M_(A,B) = mat(1,1,1,1; 1,1,1,1)$:
    - Assignment $(x_1=0, x_2=1)$ gives boundary $bold(b) = (x_1 xor x_2) = 1 in bb(F)_2^1$ (parity)
    - All 4 assignments $(0,0), (0,1), (1,0), (1,1)$ map to just 2 boundary values: $0$ or $1$
- *State aggregation*: Group by boundary equivalence, sum counts over equivalent assignments

*Advantages of recursive formulation:*
1. Natural tree structure: Each call handles one node
2. Automatic memoization: Each subtree computed once
3. Clear separation: Base cases vs. merge logic
4. Compositional: Easy to modify for different constraint types

*Complexity:* $O(n dot 2^(3r))$ where $r$ is rank width, $n$ is number of nodes. Each merge processes $<= 2^r$ states from each child, taking $O(2^(2r))$ combinations.

#theorem([Tractability via Rank Width])[
  For a CNF formula with rank width $r$, model counting can be solved in time $O(n dot 2^(3r))$, polynomial for constant $r$. For systems of XOR constraints, $r$ is bounded by the overlap degree, making such problems tractable even with dense CNF encodings.
]
