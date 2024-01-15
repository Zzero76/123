from ase.build import bulk, make_supercell
from ase.calculators.espresso import Espresso
from ase.optimize.bfgs import BFGS
from ase.constraints import ExpCellFilter

# Create a Silicon unit cell
a = 5.43  # lattice constant for Silicon in Angstrom
si = bulk('Si', 'diamond', a=a)

# Create supercell

supercell_size = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]  # Define the size of the supercell (5x5x5)
supercell = make_supercell(si, supercell_size)


pseudopotentials = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}  # Add Phosphorus pseudopotential

input_data = {
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'prefix': 'si_scf',
        'outdir': './',
    },
    'system': {
        'ecutwfc': 50,
        'ecutrho': 200,
        'occupations': 'smearing',
        'smearing': 'mv',
        'degauss': 0.01,
        'input_dft': 'PBE',
    },
    'electrons': {
        'conv_thr': 1.0e-8,
    },
}

# Set up the Quantum ESPRESSO calculator
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=('3', '3', '3'), input_data=input_data)

# Attach the calculator to the doped Silicon atoms
supercell.calc = calc

# Perform a geometry optimization
stress_constraints = ExpCellFilter(supercell)

opt = BFGS(stress_constraints, trajectory='Opt.traj', logfile='Opt.log')
opt.run(fmax=0.01)