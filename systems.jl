using DFTK
using Unitful
using UnitfulAtomic
using LinearAlgebra

function build_magnetic_moments(atoms::AbstractArray; magmoms...)
    magmoms = Dict(magmoms)
    map(atoms) do (element, positions)
        element => fill(magmoms[element.symbol], length(positions))
    end
end

function attach_pseudos(atoms::AbstractArray; pseudomap...)
    pseudomap = Dict(pseudomap)
    map(atoms) do (element, position)
        pspfile = get(pseudomap, element.symbol, nothing)
        ElementPsp(element.symbol, psp=load_psp(pspfile)) => position
    end
end

function create_supercell(lattice, atoms, supercell::Tuple)
    super = ase_atoms(lattice, atoms) * supercell
    load_lattice(super), load_atoms(super)
end
function create_supercell(lattice, atoms, supercell::Int)
    create_supercell(lattice, atoms, (1, 1, supercell))
end


"""
Iron test case (Many k-points, smallish number of bands):
    - Larger  `supercell` increases bands, increases FFT size, decreases k-points
    - Larger  `Ecut`      increases FFT size
    - Smaller `kspacing`  increases k-points
"""
function make_iron(; supercell=(1, 1, 1), Ecut=40, kspacing=0.13)
    a = 5.42352  # Bohr
    lattice = a / 2 * [[-1  1  1]; [ 1 -1  1]; [ 1  1 -1]]
    atoms = [ElementCoulomb(:Fe) => [zeros(3)]]
    lattice, atoms = create_supercell(lattice, atoms, supercell)
    atoms = attach_pseudos(atoms, Fe="hgh/pbe/fe-q16")
    magnetic_moments = build_magnetic_moments(atoms, Fe=4)

    model = model_PBE(lattice, atoms;
                      smearing=Smearing.Gaussian(),
                      temperature=1e-3,
                      magnetic_moments)
    (; basis=PlaneWaveBasis(model; Ecut, kgrid=kgrid_from_minimal_spacing(model, kspacing)),
       magnetic_moments)
end


"""
Aluminium test case (Sort of balanced):
    - Larger  `supercell` increases bands, increases FFT size, decreases k-points
    - Larger  `Ecut`      increases FFT size
    - Smaller `kspacing`  increases k-points
"""
function make_aluminium(; supercell=(3, 1, 1), Ecut=30, kspacing=0.15)
    lattice = 7.65339 * I(3)
    atoms   = [ElementCoulomb(:Al) => [[0, 0, 0],     [0, 1/2, 1/2],
                                       [1/2, 0, 1/2], [1/2, 1/2, 0]]]
    lattice, atoms = create_supercell(lattice, atoms, supercell)
    atoms = attach_pseudos(atoms, Al="hgh/pbe/al-q3.hgh")

    model = model_PBE(lattice, atoms;
                      smearing=Smearing.Gaussian(),
                      temperature=1e-3)
    (basis=PlaneWaveBasis(model; Ecut, kgrid=kgrid_from_minimal_spacing(model, kspacing)),
     magnetic_moments=[])
end


"""
Caffeine molecule test case (mostly bands and large FFT size):
    - Larger  `Ecut`     increases FFT size
    - Smaller `boxsize`  increases FFT size
"""
# 1 k-Point, large cell size, small number of bands
# (a very non-standard use case)
function make_caffeine(; Ecut=25, boxsize=20)
    @assert boxsize ≥ 15
    caffeine_xyz = """
        N    5.5015   4.8324   4.51275
        C    5.5499   3.4289   4.53495
        N    4.4279   2.7660   4.07265
        C    3.3518   3.4813   3.60165
        C    3.3514   4.8495   3.61985
        C    4.4443   5.6413   4.06195
        N    2.1172   5.1973   3.10165
        C    1.4925   4.0448   2.81165
        H    0.6190   4.0074   2.43485
        N    2.2100   2.9491   3.10105
        C    6.7257   5.5172   4.97035
        H    7.0653   5.0762   5.77575
        H    6.5188   6.4544   5.17285
        H    7.4048   5.4783   4.26445
        O    6.5268   2.8330   4.95085
        C    4.4272   1.2944   4.05585
        H    4.9878   0.9760   3.31835
        H    3.5130   0.9685   3.93275
        H    4.7825   0.9566   4.90505
        O    4.5172   6.8662   4.08205
        C    1.5560   6.5441   2.96845
        H    1.9916   7.0090   2.22425
        H    1.7073   7.0434   3.79815
        H    0.5952   6.4813   2.79415
    """
    atoms = Dict{Symbol,Vector{Vector{Float64}}}()
    for line in split(strip(caffeine_xyz), "\n")
        symbol, x, y, z = split(line)
        coords = parse.(Ref(Float64), [x, y, z]) .* austrip(1u"Å") ./ boxsize
        symbol = Symbol(symbol)
        if symbol in keys(atoms)
            push!(atoms[symbol], coords)
        else
            atoms[symbol] = [coords]
        end
    end
    atoms = [ElementCoulomb(k) => v for (k, v) in pairs(atoms)]

    lattice = boxsize * I(3)
    atoms = attach_pseudos(atoms, C="hgh/pbe/c-q4.hgh",
                                  N="hgh/pbe/n-q5.hgh",
                                  H="hgh/pbe/h-q1.hgh",
                                  O="hgh/pbe/o-q6.hgh")

    model = model_PBE(lattice, atoms)
    (basis=PlaneWaveBasis(model; Ecut, kgrid=(1, 1, 1)), magnetic_moments=[])
end
