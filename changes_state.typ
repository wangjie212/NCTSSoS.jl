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
    title: [Discussion 2025/05/27],
    subtitle: [],
    author: [Yusheng Zhao],
    date: datetime.today(),
    institution: [HKUST(GZ)],
    logo: none,
  ),
)

#set heading(numbering: numbly("{1}.", default: "1.1"))

#title-slide()

== Topic
1. Performance benchmark
2. Problems with implementation of State Polynomial Optimization

== Testcase 
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
== `cs_nctssos_first` 
#set text(14pt)
- NCTSSOS
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
- NCTSSoS.jl
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

#pause
#place(top+left, dx:120pt,dy:30pt, rect(width:400pt,fill:red.transparentize(70%)))
#place(top+left, dx:120pt,dy:220pt, rect(width:400pt,fill:red.transparentize(70%)))
#pause
#place(top+left, dx:5pt,dy:160pt, rect(width:400pt,fill:yellow.transparentize(70%)))
#place(top+left, dx:5pt,dy:350pt, rect(width:400pt,fill:yellow.transparentize(70%)))
== `cs_nctssos_higher!` 
#set text(16pt)

- NCTSSOS
```text
BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
 Single result which took 14.974 s (0.18% GC) to evaluate,
 with a memory estimate of 100.44 MiB, over 2239590 allocations.
```
- NCTSSoS.jl
 ```text
BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
 Single result which took 20.239 s (20.70% GC) to evaluate,
 with a memory estimate of 27.77 GiB, over 539154277 allocations.
 ```
#pause
#place(top+left, dx:160pt,dy:40pt, rect(height:20pt,width:190pt,fill:red.transparentize(70%)))
#place(top+left, dx:170pt,dy:140pt, rect(height:20pt,width:200pt,fill:red.transparentize(70%)))
#pause
#place(top+left, dx:200pt,dy:60pt, rect(height:20pt,width:320pt,fill:yellow.transparentize(70%)))
#place(top+left, dx:200pt,dy:160pt, rect(height:20pt,width:320pt,fill:yellow.transparentize(70%)))

==  Most Costly Part
- For each group of variables obtained from correlative sparsity, compute the term sparsity 
```julia
    cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, solver_config.ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
    end
```
- `iterate_term_sparse_supp` is the main culprit 
```text
BenchmarkTools.Trial: 77 samples with 1 evaluation per sample.
 Range (min … max):  63.343 ms … 71.434 ms  ┊ GC (min … max): 10.31% … 14.84%
 Time  (median):     64.949 ms              ┊ GC (median):    10.87%
 Time  (mean ± σ):   65.305 ms ±  1.339 ms  ┊ GC (mean ± σ):  11.03% ±  1.04%

            ▃▃ █   ▁▆       ▁      ▃           ▁               
  ▇▄▁▄▄▄▄▄▄▇██▄█▇▄▇██▇▇▁▇▄▇▄█▇▇▁▇▇▁█▄▄▄▁▁▁▄▁▁▄▁█▁▁▁▄▁▁▁▁▄▁▁▁▇ ▁
  63.3 ms         Histogram: frequency by time        68.4 ms <

 Memory estimate: 328.74 MiB, allocs estimate: 4885431.
```



== Issues
- State Polynomial Implemented and tested
- Why is basis in state polynomial the same as tracial? #link("https://github.com/wangjie212/NCTSSOS/blob/23fb8410ef4e646a6488e62bc9b03948cfce984d/src/trace.jl#L101")[code location]
- #image("statepolytracediff.png")
- See reference to "State polynomials: positivity, optimization and nonlinear Bell inequalities"


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


