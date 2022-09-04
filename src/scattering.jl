function get_model2D(materialsub, materialg, lbwidth, lbfreq, ubwidth, ubfreq, orderwidth, orderfreq, surdata_dir)
    filename = @sprintf("%s/%s %s lb %.6f %.6f ub %.6f %.6f order %d %d", surdata_dir, materialsub, materialg, lbwidth, lbfreq, ubwidth, ubfreq, orderwidth, orderfreq)
    f = open(filename,"r")
    mat = readdlm(filename, Float64)
    close(f)
    
    datareal = transpose(reshape(mat[:,5], (orderfreq+1, orderwidth+1)) )
    dataimag = transpose(reshape(mat[:,6], (orderfreq+1, orderwidth+1)) )
    
    data = complex.(datareal, dataimag)
    model2D = chebinterp(data, @SVector[lbwidth, lbfreq], @SVector[ubwidth, ubfreq])
end

function get_models1D(model2D, orderfreqPSF)
    widths = chebpoints(size(model2D.coefs,1) - 1,model2D.lb[1], model2D.ub[1]);
    freqs = reverse(chebpoints(orderfreqPSF, model2D.lb[2], model2D.ub[2]))
    models1D = [chebinterp(map(w->model2D(@SVector[w,f]), widths), model2D.lb[1], model2D.ub[1]) for f in freqs]
end

# rrule for Chebyshev polynomial functor. TODO: support chebjacobian (or explicitly don't support it)
# TODO: support x real 
function rrule(c::ChebPoly, x::AbstractVector)
    project_x = ProjectTo(x)
    y, Δy = chebgradient(c, x)
    pullback(∂y) = NoTangent(), project_x(∂y * Δy')
    y, pullback
end

