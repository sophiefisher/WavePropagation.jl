function incident_field(D, f, n, gridL, cellL)
    ω = convert(typeof(f), 2) * π * f
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + D^2)
        ℯ ^ (k * r * im) / ( convert(typeof(f),4) * π * r)
    end

    grid = range(-gridL/convert(typeof(f), 2) + convert(typeof(f), 0.5), gridL/convert(typeof(f), 2) - convert(typeof(f), 0.5), length = gridL) .* cellL
    incident = [efield(x, y) for x in grid, y in grid]
end

"""
    greens(D, f, ϵ, μ, gridL, cellL)

Returns the Green's function for
the near to far transformation.

# Arguments
- `D`: distance between near and far plane
- `f`: frequency
- `ϵ`: relative permittivity
- `μ`: relative permeability
- `gridL`: side length of E-field grid (in cells)
- `cellL`: length of each cell of grid
"""
function greens(D, f, ϵ, μ, gridL, cellL)
    ω = convert(typeof(f),2) * π * f
    n = √(ϵ*μ)
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + D^2)
        D * (convert(typeof(f),-1) + k * r * im) * ℯ ^ (k * r * im) / (convert(typeof(f),4) * π * r^3)
    end

    gridout = range(-gridL, gridL - convert(typeof(f),1), length = gridL * convert(typeof(gridL),2) ) .* cellL
    [efield(x, y) * -μ / ϵ for x in gridout, y in gridout]
end
