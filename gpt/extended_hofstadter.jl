"""
    get_h(p, q, mu, tnn, tnnn, tdnnn, kx, ky;
          peierls_sign=1,
          shell_signs=(nn=1, nnn=1, dnnn=1))

Return the `q × q` magnetic Bloch Hamiltonian of an extended square-lattice
Hofstadter model in Landau gauge `A = (0, Bx, 0)`.

The magnetic unit cell contains the orbitals

    m = 0, ..., q-1,       x_m = m + 1/2,

and the dimensionless flux per ordinary unit cell is `phi = p/q`. The hopping
shells are

    tnn     : axial nearest neighbor,          (±1, 0), (0, ±1)
    tnnn    : axial distance-two neighbor,     (±2, 0), (0, ±2)
    tdnnn   : diagonal neighbor,                (±1, ±1)

All three `t` arguments are expected to be nonnegative magnitudes. Their
intrinsic signs can be changed with `shell_signs`; for example,

    shell_signs = (nn=-1, nnn=-1, dnnn=1)

uses the conventional negative sign for the axial hoppings and a positive
diagonal hopping. By default, all matrix-element signs are positive so that
the supplied `t` values are used directly.

`mu` is inserted as the onsite matrix element. If `mu` denotes the
thermodynamic chemical potential in `H - mu*N`, pass `-mu`.

Momentum convention:

    -pi/q <= kx < pi/q,       -pi <= ky < pi.

For a straight hop `(dx,dy)` starting at sublattice `m`, the Peierls phase is

    theta(m,dx,dy) =
        peierls_sign * 2pi*(p/q)*dy*(m + 1/2 + dx/2).

Thus an upward NN hop has phase `2pi*(p/q)*(m+1/2)`, while a northeast
diagonal hop has phase `2pi*(p/q)*(m+1)`. Set `peierls_sign=-1` if the opposite
charge/Peierls convention is used.
"""
function get_h(
    p::Integer,
    q::Integer,
    mu::Real,
    tnn::Real,
    tnnn::Real,
    tdnnn::Real,
    kx::Real,
    ky::Real;
    peierls_sign::Real = 1,
    shell_signs = (nn=1, nnn=1, dnnn=1),
)
    q > 0 || throw(ArgumentError("q must be positive"))

    for (name, t) in ((:tnn, tnn), (:tnnn, tnnn), (:tdnnn, tdnnn))
        t >= 0 || throw(ArgumentError("$name must be a nonnegative magnitude"))
        isfinite(t) || throw(ArgumentError("$name must be finite"))
    end
    all(isfinite, (mu, kx, ky, peierls_sign)) ||
        throw(ArgumentError("mu, kx, ky, and peierls_sign must be finite"))

    phi = float(p) / float(q)
    H = zeros(ComplexF64, q, q)

    for m in 0:(q - 1)
        H[m + 1, m + 1] = mu
    end

    # Add one directed representative of an undirected bond and then its
    # Hermitian conjugate. This remains correct for q=1 and q=2, where
    # distinct real-space bonds can contribute to the same matrix element.
    function add_bonds!(t::Real, sign::Real, displacements)
        isfinite(sign) || throw(ArgumentError("hopping signs must be finite"))

        for (dx, dy) in displacements
            for m in 0:(q - 1)
                unwrapped_m = m + dx
                x_cell_shift = fld(unwrapped_m, q)
                m_to = mod(unwrapped_m, q)

                x_midpoint = m + 0.5 + dx / 2
                theta = peierls_sign * 2π * phi * dy * x_midpoint

                # Cell-periodic Bloch convention. A hop crossing X magnetic
                # cells and Y ordinary y-cells contributes exp[-ik·(Xq,Y)].
                bloch_phase = -(kx * q * x_cell_shift + ky * dy)
                z = sign * t * cis(theta + bloch_phase)

                H[m_to + 1, m + 1] += z
                H[m + 1, m_to + 1] += conj(z)
            end
        end
        return nothing
    end

    add_bonds!(tnn,   shell_signs.nn,    ((1, 0), (0, 1)))
    add_bonds!(tnnn,  shell_signs.nnn,   ((2, 0), (0, 2)))
    add_bonds!(tdnnn, shell_signs.dnnn,  ((1, 1), (1, -1)))

    return H
end


# Example:
#
# using LinearAlgebra
# H = get_h(
#     1, 4,             # p, q
#     0.0,              # onsite matrix element
#     1.0, 0.05, 0.12, # |t_NN|, |t_axial-distance-2|, |t_diagonal|
#     0.0, 0.0,         # kx, ky
# )
# @assert ishermitian(H)
#
# At zero flux and q=1, the positive-sign convention reduces to
#
#   H(k) = mu
#        + 2tnn   [cos(kx)  + cos(ky)]
#        + 2tnnn  [cos(2kx) + cos(2ky)]
#        + 4tdnnn cos(kx)cos(ky).
