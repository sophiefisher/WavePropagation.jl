function get_model(materialsub, materialg, lbwidth, lbfreq, ubwidth, ubfreq, orderwidth, orderfreq, surdata_dir)
    filename = @sprintf("%s/%s %s lb %.6f %.6f ub %.6f %.6f order %d %d", surdata_dir, materialsub, materialg, lbwidth, lbfreq, ubwidth, ubfreq, orderwidth, orderfreq)
    f = open(filename,"r")
    mat = readdlm(filename, Float64)
    close(f)
    datareal = transpose(reshape(mat[:,5], (orderfreq+1, orderwidth+1)) )
    dataimag = transpose(reshape(mat[:,6], (orderfreq+1, orderwidth+1)) )
    
    modelreal = chebinterp(datareal, [lbwidth, lbfreq], [ubwidth, ubfreq])
    modelimag = chebinterp(dataimag, [lbwidth, lbfreq], [ubwidth, ubfreq])
    
    surmodel(width, freq) = modelreal([width, freq]) + (im * modelimag([width, freq]) )
    
    surmodel 
    
end

# rrule for Chebyshev polynomial functor. TODO: support chebjacobian (or explicitly don't support it)
# TODO: support x real 
function rrule(c::ChebPoly, x::AbstractVector)
    project_x = ProjectTo(x)
    y, Δy = chebgradient(c, x)
    pullback(∂y) = NoTangent(), project_x(∂y * Δy')
    y, pullback
end

