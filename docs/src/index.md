# NCTSSoS

[NCTSSoS.jl](https://github.com/wangjie212/NCTSSoS) aims to provide a user-friendly and efficient tool for solving optimization problems with non-commutative/trace/state polynomials, which is based on the structured moment-SOHS hierarchy.

## Features

- Ergonomic API: Easy to use and intuitive interface for defining polynomial optimization problems!
- General Objectives: Eigenvalue, State, Tracial Polynomial Optimizations are supported!
- Correlative and Term Sparsities: Plug and play API for utilizing sparsities in reducing cost of solving optimization problems.
- Open Source: Our code is published on GitHub, feel free to contribute!

## Installation

[NCTSSoS.jl](https://github.com/wangjie212/NCTSSoS) could be installed by running

```julia
using Pkg
Pkg.add("NCTSSoS")
```

For a quick test, please refer to [Bell Inequalities Examples](@ref bell-inequalities).

To make sure everything works, you can run the tests by executing the make command in the root directory of the package.

```bash
make test
```

## Credits

The following people are involved in the development of NCTSSoS.jl:

- [Jie Wang](https://wangjie212.github.io/jiewang), Academy of Mathematics and Systems Science, Chinese Academy of Sciences.
- [Jin-Guo Liu](https://giggleliu.github.io/), Advanced Materials Thrust , The Hong Kong University of Science and Technology(Guangzhou).
- [Yusheng Zhao](https://exaclior.github.io/), Advanced Materials Thrust , The Hong Kong University of Science and Technology(Guangzhou).
- [Huanhai Zhou](https://github.com/fliingelephant), Advanced Materials Thrust , The Hong Kong University of Science and Technology(Guangzhou).

If this project is useful for your work please consider

- Citing NCTSSoS.jl
- Leave us a star on [GitHub!](https://github.com/wangjie212/NCTSSoS.jl)!

## License

NCTSSoS.jl is licensed under the MIT License.

# Related packages

- [TSSOS](https://github.com/wangjie212/TSSOS): Commutative polynomial optimization
- [ChordalGraph](https://github.com/wangjie212/ChordalGraph): Chordal graphs and chordal extentions
