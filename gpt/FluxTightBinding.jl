module FluxTightBinding

using LinearAlgebra

export TBParameters, PARAMETER_TABLE, tb_parameters, rational_flux
export tb_hamiltonian, bands, validate_model

"""
    TBParameters

Four-parameter isotropic square-lattice model in the same energy units as the
input hopping files. `tdNNN` uses the straight-line Peierls convention.
"""
struct TBParameters
    alpha::Float64
    mu::Float64
    tNN::Float64
    tNNN::Float64
    tdNNN::Float64
    origin::Symbol
    quality::Symbol
end

# Each row retains both the recommended model values and the direct extraction.
# For the failed 5/7 and 7/5 calculations, the recommended values are linear
# interpolations between the adjacent usable fluxes.
const PARAMETER_TABLE = [
    (p=0, q=1, alpha=0.0, mu=-0.031030036, tNN=0.0006277,
     tNNN=1.5045e-7, tdNNN=-1.473793319102998e-6, quality=:moderate,
     usable=true, origin=:extracted, raw_mu=-0.0309714285714,
     raw_tNN=0.000610787880, raw_tNNN=2.452837693e-5,
     raw_tdNNN=-3.243200500e-5),
    (p=2, q=7, alpha=2 / 7, mu=-0.0309714285714, tNN=0.000610787880,
     tNNN=2.452837693e-5, tdNNN=-3.243200500e-5, quality=:moderate,
     usable=true, origin=:extracted, raw_mu=-0.0309714285714,
     raw_tNN=0.000610787880, raw_tNNN=2.452837693e-5,
     raw_tdNNN=-3.243200500e-5),
    (p=3, q=7, alpha=3 / 7, mu=-0.0309114285714, tNN=0.000606190480,
     tNNN=3.609333329e-5, tdNNN=-3.798190300e-5, quality=:noisy,
     usable=true, origin=:extracted, raw_mu=-0.0309114285714,
     raw_tNN=0.000606190480, raw_tNNN=3.609333329e-5,
     raw_tdNNN=-3.798190300e-5),
    (p=1, q=2, alpha=1 / 2, mu=-0.0308075000000, tNN=0.000614718270,
     tNNN=1.683255850e-5, tdNNN=-1.711845600e-5, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0308075000000,
     raw_tNN=0.000614718270, raw_tNNN=1.683255850e-5,
     raw_tdNNN=-1.711845600e-5),
    (p=3, q=5, alpha=3 / 5, mu=-0.0307556000000, tNN=0.000577234690,
     tNNN=1.505966222e-5, tdNNN=-1.342539300e-5, quality=:moderate,
     usable=true, origin=:extracted, raw_mu=-0.0307556000000,
     raw_tNN=0.000577234690, raw_tNNN=1.505966222e-5,
     raw_tdNNN=-1.342539300e-5),
    (p=5, q=8, alpha=5 / 8, mu=-0.0307336250000, tNN=0.000586299700,
     tNNN=1.649830411e-5, tdNNN=-1.311027600e-5, quality=:usable,
     usable=true, origin=:extracted, raw_mu=-0.0307336250000,
     raw_tNN=0.000586299700, raw_tNNN=1.649830411e-5,
     raw_tdNNN=-1.311027600e-5),
    (p=5, q=7, alpha=5 / 7, mu=-0.0306355000000, tNN=0.000559786180,
     tNNN=2.328963200e-5, tdNNN=3.629079200e-6, quality=:failed,
     usable=false, origin=:safe_interpolation, raw_mu=-0.0336980000000,
     raw_tNN=0.009069297938, raw_tNNN=-0.001903105313,
     raw_tdNNN=-7.239958839e-6),
    (p=3, q=4, alpha=3 / 4, mu=-0.0305962500000, tNN=0.000549180780,
     tNNN=2.600616273e-5, tdNNN=1.032482100e-5, quality=:moderate,
     usable=true, origin=:extracted, raw_mu=-0.0305962500000,
     raw_tNN=0.000549180780, raw_tNNN=2.600616273e-5,
     raw_tdNNN=1.032482100e-5),
    (p=6, q=7, alpha=6 / 7, mu=-0.0304587142857, tNN=0.000554799740,
     tNNN=1.885414144e-5, tdNNN=-1.383209200e-5, quality=:usable,
     usable=true, origin=:extracted, raw_mu=-0.0304587142857,
     raw_tNN=0.000554799740, raw_tNNN=1.885414144e-5,
     raw_tdNNN=-1.383209200e-5),
    (p=1, q=1, alpha=1.0, mu=-0.0301940000000, tNN=0.000534510540,
     tNNN=2.177062484e-5, tdNNN=-2.498011100e-5, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0301940000000,
     raw_tNN=0.000534510540, raw_tNNN=2.177062484e-5,
     raw_tdNNN=-2.498011100e-5),
    (p=7, q=6, alpha=7 / 6, mu=-0.0299881666667, tNN=0.000493891550,
     tNNN=1.252454843e-5, tdNNN=-1.079918300e-5, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0299881666667,
     raw_tNN=0.000493891550, raw_tNNN=1.252454843e-5,
     raw_tdNNN=-1.079918300e-5),
    (p=7, q=5, alpha=7 / 5, mu=-0.0294847100000, tNN=0.000442798450,
     tNNN=1.086826600e-5, tdNNN=4.564370900e-7, quality=:failed,
     usable=false, origin=:safe_interpolation, raw_mu=-0.0299112000000,
     raw_tNN=0.000798050353, raw_tNNN=-8.935699973e-5,
     raw_tdNNN=-1.539735026e-5),
    (p=5, q=3, alpha=5 / 3, mu=-0.0289093333333, tNN=0.000384406330,
     tNNN=8.975372709e-6, tdNNN=1.332000300e-5, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0289093333333,
     raw_tNN=0.000384406330, raw_tNNN=8.975372709e-6,
     raw_tdNNN=1.332000300e-5),
    (p=7, q=4, alpha=7 / 4, mu=-0.0287000000000, tNN=0.000370079120,
     tNNN=8.292691817e-6, tdNNN=1.467292200e-5, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0287000000000,
     raw_tNN=0.000370079120, raw_tNNN=8.292691817e-6,
     raw_tdNNN=1.467292200e-5),
    (p=2, q=1, alpha=2.0, mu=-0.0279630000000, tNN=0.000319025190,
     tNNN=5.988054798e-6, tdNNN=1.299356100e-5, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0279630000000,
     raw_tNN=0.000319025190, raw_tNNN=5.988054798e-6,
     raw_tdNNN=1.299356100e-5),
    (p=8, q=3, alpha=8 / 3, mu=-0.0257720000000, tNN=0.000202925820,
     tNNN=2.999824796e-6, tdNNN=-4.188372400e-6, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0257720000000,
     raw_tNN=0.000202925820, raw_tNNN=2.999824796e-6,
     raw_tdNNN=-4.188372400e-6),
    (p=3, q=1, alpha=3.0, mu=-0.0244290000000, tNN=0.000159501810,
     tNNN=-5.229026183e-6, tdNNN=1.710368200e-6, quality=:clean,
     usable=true, origin=:extracted, raw_mu=-0.0244290000000,
     raw_tNN=0.000159501810, raw_tNNN=-5.229026183e-6,
     raw_tdNNN=1.710368200e-6),
]

const _SAFE_ROWS = filter(row -> row.usable, PARAMETER_TABLE)

function _linear_value(alpha::Real, key::Symbol; extrapolate::Bool=false)
    rows = _SAFE_ROWS
    amin, amax = rows[1].alpha, rows[end].alpha
    if !extrapolate && !(amin <= alpha <= amax)
        throw(DomainError(alpha, "flux is outside the fitted interval [$amin, $amax]"))
    end

    if alpha <= amin
        lo, hi = rows[1], rows[2]
    elseif alpha >= amax
        lo, hi = rows[end - 1], rows[end]
    else
        hi_index = findfirst(row -> row.alpha >= alpha, rows)
        hi = rows[hi_index]
        if isapprox(alpha, hi.alpha; atol=1e-14, rtol=0)
            return Float64(getproperty(hi, key))
        end
        lo = rows[hi_index - 1]
    end
    weight = (alpha - lo.alpha) / (hi.alpha - lo.alpha)
    return (1 - weight) * getproperty(lo, key) + weight * getproperty(hi, key)
end

"""
    tb_parameters(flux; mode=:safe, extrapolate=false)

Return `mu`, `tNN`, `tNNN`, and `tdNNN` at a physical flux. In `:safe` mode,
linearly interpolate only through usable source rows. Thus the failed 5/7 and
7/5 rows are automatically replaced. `mode=:raw` is available only at one of
the 16 tabulated rational fluxes and exposes the direct numerical extraction.
"""
function tb_parameters(flux::Real; mode::Symbol=:safe, extrapolate::Bool=false)
    alpha = Float64(flux)
    if mode == :raw
        index = findfirst(row -> isapprox(alpha, row.alpha; atol=1e-14, rtol=0),
                          PARAMETER_TABLE)
        isnothing(index) &&
            throw(ArgumentError("raw parameters exist only at tabulated fluxes"))
        row = PARAMETER_TABLE[index]
        return TBParameters(alpha, row.raw_mu, row.raw_tNN, row.raw_tNNN,
                            row.raw_tdNNN, :raw, row.quality)
    elseif mode != :safe
        throw(ArgumentError("mode must be :safe or :raw"))
    end

    exact = findfirst(row -> row.usable &&
                             isapprox(alpha, row.alpha; atol=1e-14, rtol=0),
                      PARAMETER_TABLE)
    origin = isnothing(exact) ? :interpolated : :extracted
    quality = isnothing(exact) ? :interpolated : PARAMETER_TABLE[exact].quality
    return TBParameters(
        alpha,
        _linear_value(alpha, :mu; extrapolate),
        _linear_value(alpha, :tNN; extrapolate),
        _linear_value(alpha, :tNNN; extrapolate),
        _linear_value(alpha, :tdNNN; extrapolate),
        origin,
        quality,
    )
end

tb_parameters(p::Integer, q::Integer; kwargs...) =
    tb_parameters(p / q; kwargs...)

"""
    rational_flux(alpha; qmax=64)

Find the closest rational `p//q` with `1 ≤ q ≤ qmax`. A finite magnetic Bloch
Hamiltonian requires rational flux.
"""
function rational_flux(alpha::Real; qmax::Integer=64)
    qmax >= 1 || throw(ArgumentError("qmax must be positive"))
    best_p, best_q = round(Int, alpha), 1
    best_error = abs(alpha - best_p)
    for q in 1:qmax
        p = round(Int, alpha * q)
        error = abs(alpha - p / q)
        if error < best_error
            best_p, best_q, best_error = p, q, error
        end
    end
    return best_p // best_q
end

function _add_directed_hop!(
    H::Matrix{ComplexF64},
    m::Int,
    dx::Int,
    dy::Int,
    amplitude::Real,
    p::Int,
    q::Int,
    Kx::Real,
    ky::Real,
)
    destination = m + dx
    cell_shift = fld(destination, q)
    mp = mod(destination, q)
    peierls = 2π * (p / q) * dy * (m + dx / 2)
    value = amplitude * cis(cell_shift * Kx + dy * ky + peierls)
    H[mp + 1, m + 1] += value
    H[m + 1, mp + 1] += conj(value)
    return H
end

"""
    tb_hamiltonian(flux::Rational, Kx, ky; mode=:safe, extrapolate=false)

Build the `q × q` magnetic Bloch Hamiltonian in Landau gauge `A=(0,Bx)`.
`Kx` is conjugate to one q-site magnetic cell and `ky` to one lattice spacing.
The model contains isotropic NN, axial NNN, and both diagonal NNN directions.

The diagonal-NNN Peierls phase is the straight-line integral,
`2π*alpha*dy*(m + dx/2)`. Changing that microscopic path convention requires
refitting `tdNNN`.
"""
function tb_hamiltonian(
    flux::Rational,
    Kx::Real,
    ky::Real;
    mode::Symbol=:safe,
    extrapolate::Bool=false,
)
    p, q = numerator(flux), denominator(flux)
    pars = tb_parameters(p / q; mode, extrapolate)
    H = zeros(ComplexF64, q, q)
    for m in 0:(q - 1)
        H[m + 1, m + 1] += pars.mu
        _add_directed_hop!(H, m, 1, 0, pars.tNN, p, q, Kx, ky)
        _add_directed_hop!(H, m, 0, 1, pars.tNN, p, q, Kx, ky)
        _add_directed_hop!(H, m, 2, 0, pars.tNNN, p, q, Kx, ky)
        _add_directed_hop!(H, m, 0, 2, pars.tNNN, p, q, Kx, ky)
        _add_directed_hop!(H, m, 1, 1, pars.tdNNN, p, q, Kx, ky)
        _add_directed_hop!(H, m, -1, 1, pars.tdNNN, p, q, Kx, ky)
    end
    return Hermitian(H)
end

tb_hamiltonian(p::Integer, q::Integer, Kx::Real, ky::Real; kwargs...) =
    tb_hamiltonian(p // q, Kx, ky; kwargs...)

"""
    tb_hamiltonian(alpha::AbstractFloat, Kx, ky; qmax=64, tol=1e-8, kwargs...)

Approximate a floating-point flux by a rational with denominator at most
`qmax`. The approximation must be within `tol`.
"""
function tb_hamiltonian(
    alpha::AbstractFloat,
    Kx::Real,
    ky::Real;
    qmax::Integer=64,
    tol::Real=1e-8,
    kwargs...,
)
    flux = rational_flux(alpha; qmax)
    abs(Float64(flux) - alpha) <= tol ||
        throw(ArgumentError("no rational flux within tol=$tol and qmax=$qmax"))
    return tb_hamiltonian(flux, Kx, ky; kwargs...)
end

"""Return the sorted magnetic subband energies."""
bands(flux, Kx::Real, ky::Real; kwargs...) =
    eigvals(tb_hamiltonian(flux, Kx, ky; kwargs...))

"""Run lightweight structural checks on all tabulated usable fluxes."""
function validate_model(; atol::Real=1e-12)
    for row in _SAFE_ROWS
        flux = row.p // row.q
        for (Kx, ky) in ((0.0, 0.0), (0.37, -1.11), (-2.2, 2.4))
            H = Matrix(tb_hamiltonian(flux, Kx, ky))
            ishermitian(H) || error("non-Hermitian matrix at flux $flux")
            maximum(abs, H - H') <= atol ||
                error("Hermiticity tolerance failed at flux $flux")
            length(eigvals(Hermitian(H))) == row.q ||
                error("wrong band count at flux $flux")
        end
    end
    return true
end

end # module
