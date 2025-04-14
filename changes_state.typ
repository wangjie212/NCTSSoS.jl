#import "@preview/touying:0.6.1": *
#import themes.university: *
#import "@preview/cetz:0.3.2"
#import "@preview/fletcher:0.5.4" as fletcher: node, edge
#import "@preview/numbly:0.1.0": numbly
#import "@preview/theorion:0.3.2": *
#import "@preview/physica:0.9.5": *
#import cosmos.clouds: *
#show: show-theorion

// cetz and fletcher bindings for touying
#let cetz-canvas = touying-reducer.with(reduce: cetz.canvas, cover: cetz.draw.hide.with(bounds: true))
#let fletcher-diagram = touying-reducer.with(reduce: fletcher.diagram, cover: fletcher.hide)

#show: university-theme.with(
  aspect-ratio: "16-9",
  // align: horizon,
  // config-common(handout: true),
  config-common(frozen-counters: (theorem-counter,)),  // freeze theorem counter for animation
  config-info(
    title: [State Polynomial Implementation],
    subtitle: [],
    author: [Yusheng Zhao],
    date: datetime.today(),
    institution: [Institution],
    logo: none,
  ),
)

#set heading(numbering: numbly("{1}.", default: "1.1"))

#title-slide()

== Outline <touying:hidden>

#components.adaptive-columns(outline(title: none, indent: 1em))

== Non-commuting State Polynomial Implementation
#fletcher-diagram(
  node-stroke: .1em,
  spacing: 10em,

	node((0,0), [`StatePolynomialOp`: \  `ncstate_terms`: #text(red)[$0.2  angle.l x y angle.r  angle.l y^2 angle.r x^2$] #text[$+ 0.1 angle.l x^2 angle.r angle.l x y angle.r x y$]], inset: 15pt),
	pause,
	// if there is no term, sorting and combining like terms is a nightmare
	edge((-0.2,0), (-0.1,0.3), "->"),
	node((-0.2,0.4), [`NCStateTerm`: \ `coef`: 0.2 \ `ncstate_word`: #text(red)[$angle.l x y angle.r angle.l y^2 angle.r x^2$ ]], inset:10pt),
	pause,
	edge((-0.1,0.55),(-0.1,0.65), "->"),
	node((0,0.8), [`NCStateWord`: \ `nc_word`: $x^2$ \ `sw`: #text(red)[$angle.l x^2 angle.r angle.l x y angle.r$ ]], inset:10pt),
	pause,
	edge((0.1,0.9),(0.3,0.9), "->"),
	node((0.5,0.9), [`StateWord`: \ `state_monos::Vector{Monomial}`: \ [$x^2 , x y$] ],inset:15pt)
)

== State Polynomial Implementation
- `expval(a::StatePolynomialOp) -> StatePolynomial`
#fletcher-diagram(
  node-stroke: .1em,
  spacing: 10em,
	node((0,0), [`StatePolynomial`: \  `state_terms`: #text(red)[$0.2  angle.l x y angle.r  angle.l y^2 angle.r angle.l x^2 angle.r $] #text[$+ 0.1 angle.l x^2 angle.r angle.l x y angle.r angle.l x y angle.r $]], inset: 15pt),
	// if there is no term, sorting and combining like terms is a nightmare
	edge((-0.2,0), (-0.1,0.3), "->"),
	node((-0.2,0.4), [`StateTerm`: \ `coef`: 0.2 \ `state_word`: #text(red)[$angle.l x y angle.r angle.l y^2 angle.r angle.l x^2 angle.r $ ]], inset:10pt),
	edge((-0.1,0.55),(0.2,0.50), "->"),
	node((0.45,0.5), [`StateWord`: \ `state_monos::Vector{Monomial}`: \ [$x y$, $y^2$, $x^2$] ], inset:10pt),
)

== Basis Problem
`NCTSSOS.jl` appear to have implemented the basis for state words as `tr`
instead of expectation value of the operator
