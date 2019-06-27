# Optimal Transport for data recoding

[Regularized optimal transport of covariates and outcomes in data recoding](https://hal.archives-ouvertes.fr/hal-02123109/file/OTRecoding.pdf), Garès, Valérie and Omer, Jérémy, 2019.

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

  Run one given method on a given number of data files of a given directory The data files must be the only files with
  extension ".txt" in the directory path: name of the directory nbfiles: number of files considered, 0 if all the data files
  are tested norme : 1 or 2, norm used for distances in the space of covariates (see runallmethods for the description of other
  parameters)
```

Copyright © 2019 Jeremy Omer <jeremy.omer@insa-rennes.fr>.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as published by
the Free Software Foundation. See LICENCE file.
