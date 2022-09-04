module WavePropagation
    
    using FFTW
    using FastChebInterp
    using FastChebInterp: ChebPoly
    using Memoize
    using ChainRulesCore: ProjectTo, NoTangent
    import ChainRulesCore.rrule
    using Printf
    using DelimitedFiles
    using StaticArrays

    export planned_fft, planned_ifft
    export convolve, convolveT, convolve!, convolveT!
    export incident_field, greens
    export get_model2D, get_models1D

    include("planned_fft.jl")
    include("fields.jl")
    include("convolve.jl")
    include("scattering.jl")

end
