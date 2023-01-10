#=
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
    (;models1D, freqs)
end
=#

function get_models1D(materialsub, materialg, in_air, lbfreq, ubfreq, orderfreq, lbwidth, ubwidth, orderwidth, surdata_dir)
    filename = @sprintf("%s/%s_%s_%s_freq_%.3f_%.3f_%d_width_%.3f_%.3f_%d",surdata_dir, materialsub, materialg, in_air, lbfreq, ubfreq, orderfreq, lbwidth, ubwidth, orderwidth)
    f = open(filename,"r")
    mat = readdlm(filename, ',',Float64)
    close(f)
    data = reshape( complex.(mat[:,1], mat[:,2]), (orderwidth+1,orderfreq+1) )
    models1D = [chebinterp(data[:,i], lbwidth, ubwidth) for i in 1:orderfreq+1]
    freqs = reverse(chebpoints(orderfreq, lbfreq, ubfreq))
    (;models1D, freqs)
end
    
function rrule(c::ChebPoly, x::Real)
    project_x = ProjectTo(x)
    y, Δy = chebgradient(c, x)
    pullback(∂y) = NoTangent(), project_x(Δy' * ∂y)
    y, pullback
end