# OTRecod.jl

## Data recoding

When databases are constructed from heterogeneous sources,
it is not unusual that different encodings are used for the same
outcome. In such case, it is necessary to recode the outcome variable
before merging two databases. The method proposed for the recoding
is an application of optimal transportation where we search for a
bijective mapping between the distributions of such variable in two
databases. In this article, we build upon the work by Garés et al.,
where they transport the distributions of categorical outcomes
assuming that they are distributed equally in the two databases.
Here, we extend the scope of the model to treat all the situations
where the covariates explain the outcomes similarly in the two
databases. In particular, we do not require that the outcomes be
distributed equally. For this, we propose a model where joint
distributions of outcomes and covariates are transported. We also
propose to enrich the model by relaxing the constraints on marginal
distributions and adding an L1 regularization term. The performances
of the models are evaluated in a simulation study, and they are
applied to a real dataset. 

Valérie Garès & Jérémy Omer (2022) Regularized Optimal Transport of Covariates and Outcomes in Data Recoding, Journal of the American Statistical Association, 117:537, 320-333, DOI: [10.1080/01621459.2020.1775615](https://doi.org/10.1080/01621459.2020.1775615)

[pdf](https://hal.archives-ouvertes.fr/hal-02123109/file/OTRecoding.pdf)


## Installation

In a Julia session switch to `pkg>` mode to add `NPSMC`:

```julia
julia>] # switch to pkg> mode
pkg> add https://github.com/otrecoding/OTRecod.jl
```

Alternatively, you can achieve the above using the `Pkg` API:

```julia
julia> using Pkg
julia> pkg"add https://github.com/otrecoding/OTRecod.jl"
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

## Instance

```@docs
Instance
```

## Solution

```@docs
Solution
```
