# OTRecod.jl

Documentation for OTRecod.jl

## Installation

In a Julia session switch to `pkg>` mode to add `NPSMC`:

```julia
julia>] # switch to pkg> mode
pkg> add https://gitlab.insa-rennes.fr/otrecoding/OTRecod.jl
```

Alternatively, you can achieve the above using the `Pkg` API:

```julia
julia> using Pkg
julia> pkg"add https://gitlab.insa-rennes.fr/otrecoding/OTRecod.jl"
```

When finished, make sure that you're back to the Julian prompt (`julia>`)
and bring `OTRecod` into scope:

```julia
julia> using OTRecod
```

You can test the package with

```julia
julia>] # switch to pkg> mode
pkg> test OTRecod
```

To run an example from a dataset

```@docs
run_directory
```

```julia
julia> using OTRecod
```
