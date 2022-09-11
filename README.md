# [ARCHIVED] 
See [Deltares/GeoArrayOps.jl](https://github.com/Deltares/GeoArrayOps.jl/) for its successor.

# GeoCost

[![Build Status](https://github.com/evetion/GeoCost.jl/workflows/CI/badge.svg)](https://github.com/evetion/GeoCost.jl/actions)

Geospatial cost (friction) operations that mimic PCRaster.
These functions should however be more Julian, extensible and scale better.

For example, the `spread` function implemented here has a threshold for maximum and scales linearly with the number of cells.

## Installation
```julia
] add GeoCost
```

## Usage
Currently there are:
- Least cost (`spread`)
- Maximum cell value in area (`areamaximum`)
