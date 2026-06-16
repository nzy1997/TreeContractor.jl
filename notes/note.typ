#import "@preview/cetz:0.4.0": canvas, draw, tree, coordinate
#import "@preview/cetz-plot:0.1.2": *
#import "@preview/ctheorems:1.1.3": *
#import "@preview/algorithmic:1.0.7"

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