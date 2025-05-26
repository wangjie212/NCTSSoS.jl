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

// == Outline <touying:hidden>

// #components.adaptive-columns(outline(title: none, indent: 1em))

== Topic
1. Performance benchmark
2. Implementation of State Polynomial Optimization

== Performance Testcase 
#set text(18pt)
```julia
    n = 5
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        f += (2x[i] + 5 * x[i]^3 + 1)^2
        f -= sum([4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] + 4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset])
        f += sum([x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset])
    end

	cons = typeof(f)[]
    for i = 1:n
        push!(cons, 1 - x[i]^2)
        push!(cons, x[i] - 1 / 3)
    end
```
== NCTSSOS `cs_nctssos_first` 
```text
BenchmarkTools.Trial: 21 samples with 1 evaluation per sample.
 Range (min … max):  237.811 ms … 264.973 ms  ┊ GC (min … max): 0.00% … 6.52%
 Time  (median):     246.529 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   247.556 ms ±   6.861 ms  ┊ GC (mean ± σ):  1.45% ± 2.31%

                    ▃▃  ▃               █                        
  ▇▇▇▁▁▇▁▁▁▁▇▁▇▁▇▁▁▇██▁▁█▁▇▁▁▁▁▁▁▁▇▁▁▁▁▁█▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  238 ms           Histogram: frequency by time          265 ms <

 Memory estimate: 31.34 MiB, allocs estimate: 746096.
```
== NCTSSOS `cs_nctssos_higher!` 
```text
BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
 Single result which took 14.974 s (0.18% GC) to evaluate,
 with a memory estimate of 100.44 MiB, over 2239590 allocations.
```

== NCTSSoS.jl `cs_nctssos` 
```text
BenchmarkTools.Trial: 2 samples with 1 evaluation per sample.
 Range (min … max):  4.868 s …   4.911 s  ┊ GC (min … max): 29.25% … 28.52%
 Time  (median):     4.890 s              ┊ GC (median):    28.88%
 Time  (mean ± σ):   4.890 s ± 30.368 ms  ┊ GC (mean ± σ):  28.88% ±  0.52%

  █                                                       █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  4.87 s         Histogram: frequency by time        4.91 s <

 Memory estimate: 8.53 GiB, allocs estimate: 128356183.
 ```

 == NCTSSoS.jl `cs_nctssos_higher`
 ```text
BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
 Single result which took 20.239 s (20.70% GC) to evaluate,
 with a memory estimate of 27.77 GiB, over 539154277 allocations.
 ```

== Issues
- State Polynomial Implemented and tested, Jie Wang's verision has bugs. Need a look at `statepolynomial.jl`

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

== Problems
- `NCTSSOS.jl` appear to have implemented the basis for state words as `tr` instead of expectation value of the operator
- Term Sparsity seems to not require squared terms, why is that? Any reference to term sparsity for State Polynomial? Squared Terms are not required in the activated basis like usual Polynomial Optimization?


