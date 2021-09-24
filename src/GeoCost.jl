module GeoCost

using FillArrays
using StaticArrays
using Distances
using OffsetArrays
using DataStructures
using Statistics

const sqrt2 = sqrt(2.0)
const neib_8 = @SMatrix[1. 1 1; 1 0 1; 1 1 1]
const distance_8 = @SMatrix[sqrt2 1 sqrt2; 1 Inf 1; sqrt2 1 sqrt2]
const distance_4 = @SMatrix[Inf 1 Inf; 1 Inf 1; Inf 1 Inf]
const Δ = CartesianIndex(1, 1)

# Write your package code here.
"""Total friction distance spread from `points`."""
function spread(points::Matrix{<:Real}, initial::Matrix{<:AbstractFloat}, friction::Matrix{<:AbstractFloat}; res=1, limit=Inf)

    ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    ofriction[begin + 1:end - 1,begin + 1:end - 1] .= friction

    result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    r = @view result[1:end - 1,1:end - 1]
    locations = points .> 0
    r[locations] .= initial[locations]

    # Construct stack for locations
    mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    mask[begin + 1:end - 1,begin + 1:end - 1] .= false

    II = CartesianIndices(size(points))
    stack = Deque{CartesianIndex}()
    for I in II[locations]
        push!(stack, I)
    end

    # Step 1: Set the distance of the starting node to 0 and the distances of all other nodes to the highest value possible.
    sdata = zeros(MMatrix{3,3})
    mcell = MMatrix{3,3}(false, false, false, false, false, false, false, false, false)

    # fcdata = zeros(MMatrix{3,3})
    # rcdata = zeros(MMatrix{3,3})

    # Step 3: For each of the active node’s adjacent neighbors, set its distance to whichever is
    # less: its current distance value or the sum of the distance of the active node plus the
    # weight of the arc from the active node to that neighbor.
    while !isempty(stack)
        spread!(stack, mask, result, ofriction, sdata, mcell, res)
    end
    r
end

function spread2(points::Matrix{<:Real}, initial::Matrix{<:AbstractFloat}, friction::Matrix{<:AbstractFloat}; res=1, limit=Inf)

    ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    ofriction[begin + 1:end - 1,begin + 1:end - 1] .= friction

    result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    r = @view result[1:end - 1,1:end - 1]
    locations = points .> 0
    r[locations] .= initial[locations]

    mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    mask[begin + 1:end - 1,begin + 1:end - 1] .= false

    minval, minidx = [0.], [CartesianIndex(1, 1)]
    x = @MMatrix zeros(3, 3)

    II = CartesianIndices(size(points))
    for I ∈ II
        patch = I - Δ:I + Δ

        rdata = view(result, patch)
        fdata = view(ofriction, patch)

        x .= (fdata .+ fdata[2,2]) .* res ./ 2 .* distance_8 .+ rdata
        findmin!(minval, minidx, x)
        rdata[2,2] = min(rdata[2,2], minval[1])
    end
    for I ∈ reverse(II)
        patch = I - Δ:I + Δ

        rdata = view(result, patch)
        fdata = view(ofriction, patch)
        x .= (fdata .+ fdata[2,2]) .* res ./ 2 .* distance_8 .+ rdata
        findmin!(minval, minidx, x)
        rdata[2,2] = min(rdata[2,2], minval[1])
    end
    for I ∈ II
        patch = I - Δ:I + Δ

        rdata = view(result, patch)
        fdata = view(ofriction, patch)

        x .= (fdata .+ fdata[2,2]) .* res ./ 2 .* distance_8 .+ rdata
        findmin!(minval, minidx, x)
        rdata[2,2] = min(rdata[2,2], minval[1])
    end
    r
end

function spread3(points::Matrix{<:Real}, initial::Matrix{<:AbstractFloat}, friction::Matrix{<:AbstractFloat}; res=1, limit=Inf)

    ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    ofriction[begin + 1:end - 1,begin + 1:end - 1] .= friction

    result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    r = @view result[1:end - 1,1:end - 1]
    locations = points .> 0
    r .= initial

    mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    mask[begin + 1:end - 1,begin + 1:end - 1] .= false

    minval, minidx = [0.], [CartesianIndex(1, 1)]
    x = @MMatrix zeros(3, 3)

    II = CartesianIndices(size(points))
    for I ∈ II
        patch = I - Δ:I + Δ

        rdata = view(result, patch)
        fdata = view(ofriction, patch)

        x .= (fdata .+ fdata[2,2]) .* res ./ 2 .* distance_8 .+ rdata
        findmin!(minval, minidx, x)
        rdata[2,2] = min(rdata[2,2], minval[1])
    end
    for I ∈ reverse(II)
        patch = I - Δ:I + Δ

        rdata = view(result, patch)
        fdata = view(ofriction, patch)
        x .= (fdata .+ fdata[2,2]) .* res ./ 2 .* distance_8 .+ rdata
        findmin!(minval, minidx, x)
        rdata[2,2] = min(rdata[2,2], minval[1])
    end
    for I ∈ II
        patch = I - Δ:I + Δ

        rdata = view(result, patch)
        fdata = view(ofriction, patch)

        x .= (fdata .+ fdata[2,2]) .* res ./ 2 .* distance_8 .+ rdata
        findmin!(minval, minidx, x)
        rdata[2,2] = min(rdata[2,2], minval[1])
    end
    r
end

function spread!(stack, mask, result, ofriction, sdata, mcell, res)
    I = popfirst!(stack)
    mask[I] = true
    patch = I - Δ:I + Δ

    rdata = view(result, patch)
    fdata = view(ofriction, patch)

    # fcdata .= fdata[2,2]
    # rcdata .= rdata[2,2]

    # New distance is cell_distance + average friction values
    for i ∈ eachindex(sdata)
        sdata[i] = muladd(fdata[i] + fdata[2,2], res / 2 * distance_8[i], rdata[2,2])
        mcell[i] = sdata[i] < rdata[i]  # cells where new distance is lower
    end
    rdata[mcell] .= sdata[mcell]
    result[patch] .= rdata

    # Add new cells to stack
    for I in patch[mcell]
        mask[I] || push!(stack, I)
    end
end

"""Optimized (and more accurate) function based on the same friction everywhere."""
function spread(points::Matrix{<:AbstractFloat}, initial::AbstractFloat, friction::AbstractFloat; distance=Euclidean(), res=1.0)
    locations = points .> 0
    I = CartesianIndices(size(points))

    result = fill(Inf, size(points))
    for location ∈ I[locations]
        for cell ∈ I
            result[cell] = min(evaluate(distance, location.I, cell.I) * res * friction + initial, result[cell])
        end
    end
    # result .+ initial
    m = .~isfinite.(points)
    result[m] = points[m]
    return result
end


"""Optimized (and more accurate) function based on the same friction everywhere."""
function spread(points::Matrix{<:AbstractFloat}, initial::Matrix{<:AbstractFloat}, friction::Real; distance=Euclidean(), res=1.0)
    locations = points .> 0
    I = CartesianIndices(size(points))

    result = fill(Inf, size(points))
    for location ∈ I[locations]
        for cell ∈ I
            result[cell] = min(evaluate(distance, location.I, cell.I) * res * friction + initial[location], result[cell])
        end
    end
    m = .~isfinite.(points)
    result[m] .= points[m]
    return result
end

"""
    roughness(dem::Matrix{<:AbstractFloat})

Roughness is the largest inter-cell difference of a central pixel and its surrounding cell, as defined in Wilson et al (2007, Marine Geodesy 30:3-35).
"""
function roughness(dem::Matrix{<:AbstractFloat})

    ex_dem = OffsetMatrix(fill(Inf, size(dem) .+ 2), UnitRange.(0, size(dem) .+ 1))
    # Update center
    ex_dem[begin + 1:end - 1,begin + 1:end - 1] .= dem
    # Set edges to mirror center
    ex_dem[begin, begin + 1:end - 1] .= dem[begin, :]
    ex_dem[end, begin + 1:end - 1] .= dem[end, :]
    ex_dem[begin + 1:end - 1, begin] .= dem[:, begin]
    ex_dem[begin + 1:end - 1, end] .= dem[:, end]
    # Set corners to mirror corners of center
    ex_dem[begin, begin] = dem[begin, begin]
    ex_dem[begin, end] = dem[begin, end]
    ex_dem[end, begin] = dem[end, begin]
    ex_dem[end, end] = dem[end, end]

    roughness = similar(dem)

    x = @MMatrix zeros(3, 3)

    @inbounds for I ∈ CartesianIndices(size(roughness))
        patch = I - Δ:I + Δ
        rdata = view(ex_dem, patch)

        x .= rdata .- rdata[2,2]
        roughness[I] = maximum(abs.(x))
    end
    roughness
end

"""
    TPI(dem::Matrix{<:AbstractFloat})

TPI stands for Topographic Position Index, which is defined as the difference between a central pixel and the mean of its surrounding cells (see Wilson et al 2007, Marine Geodesy 30:3-35).
"""
function TPI(dem::Matrix{<:AbstractFloat})

    ex_dem = OffsetMatrix(fill(Inf, size(dem) .+ 2), UnitRange.(0, size(dem) .+ 1))
    # Update center
    ex_dem[begin + 1:end - 1,begin + 1:end - 1] .= dem
    # Set edges to mirror center
    ex_dem[begin, begin + 1:end - 1] .= dem[begin, :]
    ex_dem[end, begin + 1:end - 1] .= dem[end, :]
    ex_dem[begin + 1:end - 1, begin] .= dem[:, begin]
    ex_dem[begin + 1:end - 1, end] .= dem[:, end]
    # Set corners to mirror corners of center
    ex_dem[begin, begin] = dem[begin, begin]
    ex_dem[begin, end] = dem[begin, end]
    ex_dem[end, begin] = dem[end, begin]
    ex_dem[end, end] = dem[end, end]

    tpi = similar(dem)

    x = @MMatrix zeros(3, 3)

    @inbounds for I ∈ CartesianIndices(size(tpi))
        patch = I - Δ:I + Δ
        rdata = view(ex_dem, patch)

        x .= rdata .* neib_8
        tpi[I] = rdata[2,2] - mean(x)
    end
    tpi
end

"""
    TRI(dem::Matrix{<:AbstractFloat})

TRI stands for Terrain Ruggedness Index, which measures the difference between a central pixel and its surrounding cells.
This algorithm uses the square root of the sum of the square of the difference between a central pixel and its surrounding cells.
This is recommended for terrestrial use cases.
"""
function TRI(dem::Matrix{<:AbstractFloat})

    ex_dem = OffsetMatrix(fill(Inf, size(dem) .+ 2), UnitRange.(0, size(dem) .+ 1))

    # Update center
    ex_dem[begin + 1:end - 1,begin + 1:end - 1] .= dem

    # Set edges to mirror center
    ex_dem[begin, begin + 1:end - 1] .= dem[begin, :]
    ex_dem[end, begin + 1:end - 1] .= dem[end, :]
    ex_dem[begin + 1:end - 1, begin] .= dem[:, begin]
    ex_dem[begin + 1:end - 1, end] .= dem[:, end]

    # Set corners to mirror corners of center
    ex_dem[begin, begin] = dem[begin, begin]
    ex_dem[begin, end] = dem[begin, end]
    ex_dem[end, begin] = dem[end, begin]
    ex_dem[end, end] = dem[end, end]

    tri = similar(dem)

    x = @MMatrix zeros(3, 3)

    @inbounds for I ∈ CartesianIndices(size(tri))
        patch = I - Δ:I + Δ
        rdata = view(ex_dem, patch)

        x .= (rdata .- rdata[2,2]).^2
        tri[I] = sqrt(sum(x))
    end
    tri
end

"""
```
B, flags = pmf(A; ωₘ, slope, dhₘ, dh₀, cellsize)
```
Applies the progressive morphological filter by [Zhang et al. (2003)] to `A`.
# Output
- `B::Array{T,2}` Maximum allowable values
- `flags::Array{Float64,2}` A sized array with window sizes if filtered, zero if not filtered.
Afterwards, one can retrieve the resulting mask for `A` by `A .<= B` or `flags .== 0.`.
# Arguments
- `A::Array{T,2}` Input Array
- `ωₘ::Float64=20.` Maximum window size [m]
- `slope::Float64=0.01` Terrain slope [m/m]
- `dhₘ::Float64=2.5` Maximum elevation threshold [m]
- `dh₀::Float64=0.2` Initial elevation threshold [m]
- `cellsize::Float64=1.` Cellsize in [m]
[Zhang et al. (2003)] Zhang, Keqi, Shu-Ching Chen, Dean Whitman, Mei-Ling Shyu, Jianhua Yan, and Chengcui Zhang. “A Progressive Morphological Filter for Removing Nonground Measurements from Airborne LIDAR Data.” IEEE Transactions on Geoscience and Remote Sensing 41, no. 4 (2003): 872–82. [https://doi.org/10.1109/TGRS.2003.810682].
"""
function pmf(A::Array{T,2};
    ωₘ::Float64=20.,
    slope::Float64=0.01,
    dhₘ::Float64=2.5,
    dh₀::Float64=0.2,
    cellsize::Float64=1.0) where T <: Real

    # Compute windowsizes and thresholds
    ωₘ = round(Int, ωₘ / cellsize)
    κ_max = floor(Int, log2(ωₘ - 1))  # determine # iterations based on exp growth
    windowsizes = Int.(exp2.(1:κ_max)) .+ 1

    # Compute tresholds
    dwindows = vcat(windowsizes[1], windowsizes)  # prepend first element so we get 0 as diff
    window_diffs = [dwindows[i] - dwindows[i - 1] for i in 2:length(dwindows)]
    height_tresholds = [min(dhₘ, slope * window_diff * cellsize + dh₀) for window_diff in window_diffs]

    # Set up arrays
    Af = copy(A)  # array to be morphed
    nan_mask = isnan.(Af)
    Af[nan_mask] .= Inf  # Replace NaN with Inf, as to always filter these

    B = copy(A)  # max_elevation raster
    out = copy(A)  # max_elevation raster

    flags = zeros(size(A))  # 0 = ground, other values indicate window size
    flags[nan_mask] .= NaN

    mask = falses(size(A))

    # Iterate over window sizes and height tresholds
    for (ωₖ, dhₜ) in zip(windowsizes, height_tresholds)
        opening!(Af, ωₖ, out)
        for I in eachindex(A)
            mask[I] = (A[I] - Af[I]) > dhₜ
        end
        for I in eachindex(flags)
            if mask[I] && flags[I] == 0
                flags[I] = ωₖ
            end
        end
        # @info "PMF with window size $ωₖ and threshold $dhₜ filters $(sum(mask)) cells."
        B .= min.(B, Af .+ dhₜ)
    end

    B, flags
end

# First discussed here https://github.com/JuliaImages/ImageFiltering.jl/issues/179
function mapwindow!(f, img, window, out)
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(ntuple(x -> window ÷ 2, ndims(img)))
    @inbounds @simd for I in R
        patch = max(I_first, I - Δ):min(I_last, I + Δ)
        out[I] = f(view(img, patch))
    end
    out
end

"""Apply the opening operation to `A` with window size `ω`."""
function opening!(A::Array{T,2}, ω::Integer, out::Array{T,2}) where T <: Real
    mapwindow!(minimum, A, ω, out)  # erosion
    mapwindow!(maximum, out, ω, A)  # dilation
    A
end


export spread, spread2, spread3, roughness, TPI, TRI, pmf

end  # module
