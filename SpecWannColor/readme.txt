This is a set of files, written in julia, which calculate the spectrum of Landau levels in a cosine periodic potential in 2D,
calculate its topological properties based on a generated Wannier plot, and color the gaps in the spectrum according to Chern.

How to run:
typical run command from terminal
julia main_SCWC.jl "[0.05, 1.1, 0.015, 50, 499, 3, 0.4]"
where the arguments in the vector are
1. starting flux (x axis lower limit)
2. final flux (x axis upper limit)
3. strength of the potential U0 in eV
4. periodicity of the potential (unit length) in Ansgstrom
5. p - resolution of the calculation
6. NLL - n of the last Landau level included
7. factor for gap recognition (check code for more details)

Note: time for calculation is approx propto p^3, NLL^2


For running on a cluster:
Execute sub_main_SCWC.jl in frontend node, which generates and submits job files to computational nodes.
