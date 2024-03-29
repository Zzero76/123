from ase.calculators.espresso import Espresso
from ase.dft import DOS
import matplotlib.pyplot as plt
from ase.io import read


# Load structure
si = read('Opt.traj@-1')
#Identify the index of a Silicon atom to replace with Phosphorus
index_to_replace = 0  # Change this index to the desired Silicon atom index

# Create a copy of the Silicon unit cell
si_doped = si.copy()

# Replace the Silicon atom at the specified index with a Phosphorus atom
si_doped[index_to_replace].symbol = 'B'


# Need to add the new pseudopotential

pseudopotentials = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF', 'B': 'b_pbe_v1.4.uspp.F.UPF'}


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
        'smearing': 'gaussian',
        'degauss': 0.01,
        'input_dft': 'PBE',
    },
    'electrons': {
        'conv_thr': 1.0e-8,
    },
}

# Set up the Quantum ESPRESSO calculator
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(12, 12, 12), input_data=input_data)

# Attach the calculator to the Silicon atoms
si_doped.calc = calc

# Perform SCF calculation
si_doped.get_potential_energy()

# Calculate DOS
dos = DOS(calc, width=0.2)
dos.get_dos()

# Get energies and DOS values
energies = dos.get_energies()
dos_values = dos.get_dos()


# Plot DOS
plt.figure(figsize=(8, 6))
plt.plot(energies, dos_values, color='red')
plt.xlabel('Energy (eV)')
plt.ylabel('Density of States')
plt.title('Electronic DOS of Silicon')

# Plot Fermi level line at 0 energy
fermi_level = calc.get_fermi_level()
# For testing
#print(fermi_level)
plt.axvline(x=0, color='k', linestyle='--', label=f'Fermi Level: {fermi_level:.2f} eV')


plt.legend()
plt.savefig('DOS.png')
