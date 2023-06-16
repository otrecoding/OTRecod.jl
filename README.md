# OTRecod.jl

[![CI](https://github.com/otrecoding/OTRecod.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/otrecoding/OTRecod.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/otrecoding/OTRecod.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/otrecoding/OTRecod.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://otrecoding.github.io/OTRecod.jl/dev)

[Regularized optimal transport of covariates and outcomes in data recoding](https://hal.archives-ouvertes.fr/hal-02123109/file/OTRecoding.pdf), Garès, Valérie and Omer, Jérémy, 2019.

Valérie Garès & Jérémy Omer, 2022.  "Regularized Optimal Transport of Covariates and Outcomes in Data Recoding," Journal of the American Statistical Association, Taylor & Francis Journals, vol. 117(537), pages 320-333, January.

Abstract:  When databases are constructed from heterogeneous sources,
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

Keywords: 
<https://ideas.repec.org/a/taf/jnlasa/v117y2022i537p320-333.html>

## Installation

The package runs on julia 1.1 and above.
In a Julia session switch to `pkg>` mode to add the package:

```julia
julia>] # switch to pkg> mode
pkg> add https://github.com/otrecoding/OTRecod.jl
```

Alternatively, you can achieve the above using the `Pkg` API:

```julia
julia> import Pkg
julia> Pkg.add(url = "https://github.com/otrecoding/OTRecod.jl")
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

```julia
julia> using OTRecod

help?> run_directory
search: run_directory

  run_directory(path, method; outname="result.out",
                              maxrelax=0.0,
                              lambda_reg=0.0,
                              nbfiles=0,
                              norme=0,
                              percent_closest=0.2)

  Run one given method on a given number of data files of a given directory. The data files must be the only files with
  extension ".txt" in the directory.

 - `path`   : name of the directory
 - `method` : `:group` or `:joint`
 - `maxrelax`: maximum percentage of deviation from expected probability masses
 - `lambda_reg`: coefficient measuring the importance of the regularization term
 - `nbfiles`: number of files considered, 0 if all the data files are tested
 - `norme`  : 0, 1 or 2, norm used for distances in the space of covariates
 - `percent_closest`: percent of closest neighbors taken in the computation of the costs (both distance and regularization related)
 - `observed`: if nonempty, list of indices of the observed covariates; this allows to exclude some latent variables.
```

Copyright © 2020 Jeremy Omer <jeremy.omer@insa-rennes.fr>.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as published by
the Free Software Foundation. See [LICENSE file](LICENSE).
