# CLSolvers.jl

This is the code generating the simulations done in: https://arxiv.org/abs/2105.02735

To run the code you need a version of Julia installed, then you can make separate scripts or
follow the Main.jl file which can be run line by line using the Julia vscode extension.

## Instantiate

To initialize the project run these comments inside the Julia REPL (From inside the project directory)
'''
    import Pkg
    Pkg.activate(".")
    Pkg.instantiate()
'''
For more information see: https://docs.julialang.org/en/v1/stdlib/Pkg/

Now all dependencies should be downloaded and the code is ready to be run.
