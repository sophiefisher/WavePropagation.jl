function incident_field(D, f, n, gridL, cellL)
    ω = 2π * f
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + D^2)
        ℯ ^ (k * r * im) / (4π * r)
    end

    grid = range(-gridL/2 + 1/2, gridL/2 - 1/2, length = gridL) .* cellL
    [efield(x, y) for x in grid, y in grid]
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
    ω = 2π * f
    n = √(ϵ*μ)
    k = n * ω

    function efield(x, y)
        r = √(x^2 + y^2 + D^2)
        D * (-1 + k * r * im) * ℯ ^ (k * r * im) / (4 * π * r^3)
    end

    gridout = range(-gridL, gridL - 1, length = gridL * 2) .* cellL
    g = [efield(x, y) * -μ / ϵ for x in gridout, y in gridout]

    g
end
