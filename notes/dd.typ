= Decision Diagrams for SAT Solving

== Introduction

Decision diagrams are powerful data structures for representing and manipulating Boolean functions. In the context of SAT (Boolean satisfiability) solving, decision diagrams provide an alternative approach to traditional DPLL-based and CDCL (Conflict-Driven Clause Learning) methods by explicitly representing the solution space of Boolean formulas.

== Binary Decision Diagrams (BDDs)

=== Definition

A Binary Decision Diagram (BDD) is a directed acyclic graph (DAG) that represents a Boolean function $f: {0,1}^n -> {0,1}$. Each non-terminal node is labeled with a Boolean variable and has two outgoing edges (low/0-edge and high/1-edge), while terminal nodes are labeled with Boolean constants (0 or 1).

*Ordered BDDs (OBDDs):* Variables appear in a consistent order along every path from root to terminal.

*Reduced Ordered BDDs (ROBDDs):* A canonical form where:
- No two nodes have identical sub-diagrams
- No node has both edges pointing to the same target

=== Properties

- *Canonicity*: For a given variable ordering, the ROBDD representation of a function is unique
- *Compactness*: Many practical Boolean functions have polynomial-size BDDs
- *Efficient Operations*: Apply operation (AND, OR, XOR, etc.) can be performed in time $O(|B_1| times |B_2|)$ where $|B_i|$ is the size of BDD $i$

=== SAT Solving with BDDs

To solve SAT using BDDs:

1. *Construction*: Build a BDD representation of the CNF formula by converting each clause to a BDD and combining them with AND operations

2. *Satisfiability Check*: A formula is satisfiable if and only if its BDD is not the constant 0 function

3. *Solution Extraction*: Any path from root to the 1-terminal represents a satisfying assignment

*Advantages:*
- Can efficiently handle certain structured formulas
- Provides compact representation for many practical problems
- Natural support for model counting and enumeration

*Limitations:*
- Variable ordering is crucial - poor orderings can cause exponential blowup
- BDD size can grow exponentially for some formula families
- Construction phase can be expensive for large formulas

== Zero-Suppressed Decision Diagrams (ZDDs)

=== Definition

ZDDs are a variant of BDDs optimized for representing sparse sets of combinations. The key difference is the reduction rule: a node is eliminated if its high-edge (1-edge) points to the 0-terminal.

=== Applications to SAT

ZDDs are particularly useful for:
- *Solution Enumeration*: Representing all satisfying assignments compactly
- *Sparse Formulas*: When most variables don't appear in the solution
- *Constraint Programming*: Representing feasible solution sets

=== Comparison with BDDs

| Aspect | BDD | ZDD |
|--------|-----|-----|
| Reduction Rule | Both edges equal | High edge to 0 |
| Best For | Dense functions | Sparse sets |
| SAT Applications | General formulas | Solution enumeration |

== d-DNNFs (Deterministic Decomposable Negation Normal Form)

=== Structure

A more general decision diagram format where:
- Formula is in Negation Normal Form (negations only on literals)
- *Decomposability*: AND nodes have disjoint variable sets in different children
- *Determinism*: OR nodes have mutually exclusive children

=== Advantages for SAT

- More expressive than BDDs - can represent some functions exponentially more compactly
- Efficient model counting: count can be computed in linear time
- Supports tractable probabilistic inference
- Knowledge compilation: compile CNF to d-DNNF offline, then answer queries efficiently

=== Compilation Approaches

1. *Top-down*: Use DPLL-style search with component caching
2. *Bottom-up*: Start from clauses and combine using AND/OR gates
3. *Hybrid*: Combine different techniques based on structure

== Algebraic Decision Diagrams (ADDs)

ADDs extend BDDs to represent functions with non-Boolean ranges: $f: {0,1}^n -> cal(R)$ where $cal(R)$ is typically real numbers or integers.

*SAT Applications:*
- Weighted SAT (MaxSAT, weighted model counting)
- Probabilistic reasoning
- Optimization problems

== Sentential Decision Diagrams (SDDs)

=== Definition

SDDs are a more structured form of d-DNNFs that respect a given vtree (variable tree) structure. Each decision node represents a disjunction of prime-sub pairs $(p_i, s_i)$ where:
- Primes $p_i$ are mutually exclusive and exhaustive
- Each sub $s_i$ uses variables disjoint from the corresponding prime

=== Properties

- *Polytime Operations*: Apply, condition, and quantify operations are polynomial in SDD size
- *Succinctness*: Can be exponentially more succinct than BDDs
- *Canonicity*: Unique representation for given vtree
- *Tractable*: Supports polynomial-time queries

=== SAT Solving with SDDs

Similar to d-DNNFs but with better structural guarantees:
1. Choose or learn a good vtree
2. Compile CNF formula to SDD
3. Query the compiled representation

== Practical Considerations

=== Variable Ordering

For BDD-based approaches, variable ordering is critical:

*Heuristics:*
- *FORCE*: Minimize expected span
- *Linear arrangement*: Based on formula structure
- *Dynamic reordering*: Adjust ordering during construction
- *Hypergraph partitioning*: Minimize edge cuts

=== When to Use Decision Diagrams for SAT

*Favorable Scenarios:*
- Formulas with special structure (e.g., planning, verification)
- Need for multiple queries on same formula
- Model counting or enumeration required
- Probabilistic inference

*Unfavorable Scenarios:*
- Random or unstructured formulas
- Very large industrial instances
- Single satisfiability query

=== Hybrid Approaches

Modern SAT solvers sometimes combine techniques:
- Use CDCL for hard subproblems
- Build decision diagrams for structured components
- Cache learned information in diagram form

== Comparison with Traditional SAT Solvers

| Aspect | Decision Diagrams | CDCL Solvers |
|--------|-------------------|--------------|
| Representation | Explicit (diagram) | Implicit (search) |
| Memory | Can be large | Typically smaller |
| Repeated Queries | Very fast | Must re-solve |
| Model Counting | Efficient | Requires extensions |
| Scalability | Depends on structure | Generally better |
| Best Use | Structured problems | General instances |

== Applications Beyond SAT

Decision diagrams compiled from SAT formulas enable:
- *Formal Verification*: Hardware and software model checking
- *Configuration*: Product line configuration spaces
- *Planning*: AI planning with state space representation
- *Network Reliability*: Analyzing fault tolerance
- *Bioinformatics*: Analyzing genetic pathways

== Tools and Libraries

*BDD Libraries:*
- CUDD (Colorado University Decision Diagram)
- BuDDy
- JavaBDD / JDD
- Sylvan (parallel BDDs)

*Knowledge Compilation:*
- c2d (compiles CNF to d-DNNF)
- Dsharp (d-DNNF compiler)
- SDD Package (UCLA)

*Integration with SAT:*
- ABC (synthesis and verification)
- NuSMV (model checking)

== Conclusion

Decision diagrams provide a complementary approach to traditional SAT solving, particularly valuable when:
1. Formula has exploitable structure
2. Multiple queries are needed
3. Solution enumeration or counting is required
4. Integration with probabilistic reasoning is desired

The choice between decision diagram approaches and traditional SAT solvers depends on the problem structure, query patterns, and performance requirements. Hybrid approaches that combine both paradigms are an active area of research.



