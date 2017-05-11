__precompile__()

module SpatiallyFilteredViews

using ImageFiltering, ImageAxes

using Base: ViewIndex
using ImageFiltering: ReshapedIIR

export SpatiallyFilteredView

immutable SpatiallyFilteredView{T,N,P<:AxisArray,A,F} <: StreamingTimeArray{T,N}
    filter!::F
    slice::A
    parent::P
end
function (::Type{SpatiallyFilteredView{T}}){T,N,_}(filter!::Function, parent::AxisArray{_,N})
    spax = space_axes(parent)
    spinds = space_indices(parent)
    slice = AxisArray(ReadOnlyArray(similar(parent.data, T, spinds)), spax)
    SpatiallyFilteredView{T,N,typeof(parent),typeof(slice),typeof(filter!)}(filter!, slice, parent)
end

Base.indices(S::SpatiallyFilteredView) = indices(S.parent)
Base.size(S::SpatiallyFilteredView)    = size(S.parent)

AxisArrays.axes(S::SpatiallyFilteredView)         = axes(S.parent)

function Base.view(S::SpatiallyFilteredView, tindex::Axis{:time,Int})
    Pv = view(S.parent, tindex)
    S.filter!(S.slice.data.parent, Pv)
    S.slice
end
function Base.view{R<:Range}(S::SpatiallyFilteredView, tindex::Axis{:time,R})
    SpatiallyFilteredView{eltype(S)}(S.filter!, view(S.parent, tindex))
end
function Base.view{T,N}(S::SpatiallyFilteredView{T,N}, I::Vararg{ViewIndex,N})
    tind = ImageAxes.filter_time_axis(axes(S.parent), I)
    sind = ImageAxes.filter_space_axes(axes(S.parent), I)
    tslice = timeaxis(S)(tind[1])
    view(view(S, tslice), sind...)
end

space_axes(parent) = ImageAxes.filter_space_axes(axes(parent), axes(parent))
space_indices(parent) = ImageAxes.filter_space_axes(axes(parent), indices(parent))

# Filters
function highpass_trunc!{N}(dest, src, kernel::NTuple{N,ReshapedIIR})
    imfilter!(dest, src, kernel, NA())
    z = zero(eltype(dest))
    for I in eachindex(dest, src)
        dest[I] = max(z, src[I]-dest[I])
    end
    dest
end

end # module
