# Recoding data using optimal transport

## Installation


In a Julia session switch to `pkg>` mode to add `NPSMC`:

```julia
julia>] # switch to pkg> mode
pkg> add https://gitlab.insa-rennes.fr/otrecoding/OTRecoding.jl
```

Alternatively, you can achieve the above using the `Pkg` API:

```julia
julia> using Pkg
julia> pkg"add https://gitlab.insa-rennes.fr/otrecoding/OTRecoding.jl"
```

When finished, make sure that you're back to the Julian prompt (`julia>`)
and bring `OTRecoding` into scope:

```julia
julia> using OTRecoding
```

You can test the package with

```julia
julia>] # switch to pkg> mode
pkg> test OTRecoding
```

Copyright © 2019 Jeremy Omer <jeremy.omer@insa-rennes.fr>.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as published by
the Free Software Foundation. See LICENCE file.
