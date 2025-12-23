# landau-butterflies
*Writen in julia 1.12.1+0.x64.linux.gnu*

Code in julia that considers electrons in a cosine periodic potential in a large magnetic field.

## Running instructions
Working files to be executed are in the *mains/* and *short/* folder.
If you want to submit to cluster, use the submission files from the */subs/* folder instead.

When running, the packages need to be activated before being used. So execute the code from Terminal like 
>   julia --project=/path/to/MyProject main_X.jl

which tells julia to look for packages in the .toml files.

### running on the cluster
make sure to initialize packages from the .toml files on the cluster version by running in julia on the cluster

>   julia> cd("path/to/project")
>
>   julia> ]
>
>   pkg> activate .
>
>   pkg> resolve
>
>   pkg> instantiate

## main_S.jl
This is the main file, which can be executed to compute the energy spectrum and optionally perform a Wannier analysis and color the gaps in the spectrum according to Chern number.**EDIT THE OUTPUT FOLDER DIRECTLY IN THE FILE!**
It can take multiple options as arguments. These are all described by calling a help function with 
>   julia --project=/path/to/MyProject main_X.jl -h

The mandatory positional arguments are:
1. U0 - potential strength in [eV]
2. a - lattice constant of potential [nm]
3. LLmax - maximum index of Landau levels used in calculation (complexity scales as O(LLmax^3)above 100 it becomes extremely heavy; needed only at very small fluxes)
4. q - sets a maximum denominator of the flux=p/q (simplified when possible; sets resolution of the spectrum; O(q^4))
5. phi_s - starting flux
6. phi_f - final flux

The optional arguments are:
* --XLL k - restricts output to energies, belonging to the **k** lowest Landau levels 
* --XBF k - restricts output to energies, belonging to the **k** lowest Hofstadter butterflies; if none are specified, output full spectrum
* --plot (-p) - output as a plot
* --data (-d) - output as a .npz data file
* --wannier (-w) - perform Wannier analysis; if plotting, colors the spectrum according to Chern numbers of gaps;
* --plotW (-z) - also output as a plot (if -p) a Wannier diagram of the spectrum
* --varyNLL (-l) - at each flux, vary the number of Landau levels used to save computational power at higher fluxes, where less LLs are needed

### Example usage
>   julia --project=/path/to/MyProject main_S.jl 0.02 5 50 240 0.25 2.0 --XBF 1 -p -w -l -z


