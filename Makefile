JL = julia --project

default: init test

init:
	$(JL) -e 'using Pkg; Pkg.precompile()'

init-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.develop(path="."); Pkg.precompile()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile()'

update-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.update(); Pkg.precompile()'

test:
	$(JL) -e 'using Pkg; Pkg.test()'

servedocs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs = ["docs/src/assets", "docs/src/generated"])'

clean:
	rm -rf docs/build
	find . -name "*.cov" -type f -print0 | xargs -0 /bin/rm -f

test-FastPoly:
	$(JL) -e 'using Pkg; Pkg.status(); include("test/fastpoly_test/runtests.jl")'

init-bench:
	$(JL) -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("AirspeedVelocity"); Pkg.build("AirspeedVelocity")'

bench:
	benchpkg \
	--script "benchmark/benchmarks.jl"

bench-plot:
	benchpkgplot  \
    --format=pdf \
    --npart=5

.PHONY: init test
