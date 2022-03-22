@memoize function make_plan(size::Tuple)
    plan_fft(zeros(ComplexF64, size), flags=FFTW.MEASURE)
end                                            

@memoize function make_plan!(size::Tuple)
    plan_fft!(zeros(ComplexF64, size), flags=FFTW.MEASURE)
end                                            

function planned_fft(x)
    plan = make_plan(size(x))
    plan * x
end

function planned_ifft(x)
    plan = make_plan(size(x))
    plan \ x
end

function planned_fft!(x)
    plan = make_plan!(size(x))
    plan * x
end

function planned_ifft!(x)
    plan = make_plan!(size(x))
    plan \ x
end

function rrule(::typeof(planned_fft), x)
    N = prod(size(x))
    function pullback(∂x)
        NoTangent(), N * planned_ifft(∂x)
    end
    return planned_fft(x), pullback
end

function rrule(::typeof(planned_ifft), x)
    N = prod(size(x))
    function pullback(∂x)
        NoTangent(), planned_fft(∂x) / N
    end
    return planned_ifft(x), pullback
end
