"""
    convolve(inp, kernel)

Convolves an inpL x inpL array with the FFT of a centered kernel 
of size kerL x kerL to produce an output of size (kerL - inpL) x (kerL - inpL).

TODO: maybe this should return a real array
"""
function convolve(inp, kernel, plan)
    #plan must be for the size of the kernel
    inpL = size(inp)[1]
    kerL = size(kernel)[1]
    outL = kerL - inpL

    arr_pad = [inp zeros(inpL, outL); zeros(outL, inpL) zeros(outL, outL)]
    #out_pad = P \ ((P * arr_pad) .* kernel)
    out_pad = plan \ (  (plan * arr_pad) .* kernel)
    out = out_pad[inpL+1:kerL, inpL+1:kerL]

    out
end

function convolveT(out, kernel, plan)
    outL = size(out)[1]
    kerL = size(kernel)[1]
    inpL = kerL - outL

    out_pad = [zeros(inpL, inpL) zeros(inpL, outL); zeros(outL, inpL) out]
    #arr_pad = P * ((P \ out_pad) .* kernel)
    arr_pad = plan * ((plan \ out_pad) .* kernel )
    arr = arr_pad[1:inpL, 1:inpL]

    arr
end

#=
function convolveT(out, kernel)
    outL = size(out)[1]
    kerL = size(kernel)[1]
    inpL = kerL - outL

    out_pad = [zeros(inpL, inpL) zeros(inpL, outL); zeros(outL, inpL) out]
    #arr_pad = P * ((P \ out_pad) .* kernel)
    arr_pad = planned_fft(planned_ifft(out_pad) .* kernel)
    arr = arr_pad[1:inpL, 1:inpL]

    arr
end

# padded stores the intermediate padded result
function convolve!(inp, kernel, padded)
    inpL = size(inp)[1]
    kerL = size(kernel)[1]
    outL = kerL - inpL

    padded .= 0
    padded[1:inpL, 1:inpL] .= inp
    planned_fft!(padded)
    padded .*= kernel
    planned_ifft!(padded)
    out = @view padded[inpL+1:kerL, inpL+1:kerL]

    out
end

function convolveT!(out, kernel, padded)
    outL = size(out)[1]
    kerL = size(kernel)[1]
    inpL = kerL - outL

    padded .= 0
    padded[inpL+1:inpL+outL, inpL+1:inpL+outL] .= out
    planned_ifft!(padded)
    padded .*= kernel
    planned_fft!(padded)
    inp = @view padded[1:inpL, 1:inpL]

    inp
end
    
=#
          
