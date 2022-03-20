# LinearSolvers

[![Build Status](https://github.com/subinbg/LinearSolvers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/subinbg/LinearSolvers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/subinbg/LinearSolvers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/subinbg/LinearSolvers.jl)

A Julia package to iteratively solve linear problems with custom matrix-vector operations.   
A user is expected to define the following functions:
- `LinearSolvers.mul!`
- `LinearSolvers.sub!`
- `LinearSolvers.add!`
- `LinearSolvers.dot`
- `LinearSolvers.norm`
- `LinearSolvers.add_scmul!`
- `LinearSolvers.sub_scmul!`
- `LinearSolvers.eltype`
- `LinearSolvers.size`
- `LinearSolvers.fill!`
- `LinearSolvers.rand!`
- `LinearSolvers.copyto!`
`src/LinearSolvers.jl` includes examples of how to properly define those functions.  
A user may also check `test/runtests.jl`.