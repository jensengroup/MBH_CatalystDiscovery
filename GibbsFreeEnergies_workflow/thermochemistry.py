import numpy as np
from tabulate import tabulate

# CONVERSION FACTORS

JOULES2CAL = 0.2390057361
HARTREE2JOULES = 4.3597482e-18
HARTEE2CALMOL = 627509.5
AMU2KG = 1.6605402e-27
ATM2PASCAL = 101325

# CONSTANTS

GAS_CONSTANT = 8.31446261815324  # J/(K x mol)
BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
PLANCK_CONSTANT = 6.626070150e-34  # J*s
SPEED_OF_LIGHT = 299792458  # m/s
AVOGARDOS_CONSTANT = 6.02214076e23  # 1/mol


def get_indices(lines, pattern, stop_pattern=None):
    indices = []
    for i, l in enumerate(lines):
        if pattern in l:
            indices.append(i)
        if stop_pattern and stop_pattern in l:
            break
    return indices


def get_frequencies(lines):
    indices = get_indices(lines, pattern="Frequencies --")
    frequencies = []
    for idx in indices:
        line = lines[idx]
        frequencies.extend([float(f) for f in line.split()[2:]])
    return np.array(frequencies)


def get_electronic_energy(lines):
    indices = get_indices(lines, pattern="SCF Done:")
    line = lines[indices[-1]]
    energy = float(line.split()[4])
    return energy


def get_solv_energy(lines):
    try:
        indices = get_indices(lines, pattern="SMD-CDS")
        line = lines[indices[-1]]
        energy = float(line.split()[-1])
    except:
        energy = np.nan
    return energy


def get_rotational_entropy(lines):
    indices = get_indices(lines, pattern="Rotational")
    line = lines[indices[-2]]
    Sr = float(line.split()[-1])
    return Sr


def get_temperature_and_pressure(lines):
    indices = get_indices(lines, "Temperature")
    line = lines[indices[0]]
    T, p = line.split()[1::3]
    return float(T), float(p)


def get_molar_mass(lines):
    indices = get_indices(lines, pattern="Molecular mass:")
    line = lines[indices[-1]]
    M = float(line.split()[2])
    return M


def clip_frequencies(frequencies, f_cutoff, verbose=True):
    """Clips Frequencies below a cutoff value to NaN
    Parameters:
    frequencies : array
        Array/List of frequencies in cm^-1
    f_cutoff : float
        Frequency Cutoff in cm^-1
    Returns:
        Clipped Frequencies
    """

    frequencies = np.asarray(frequencies)
    imag_freq = frequencies <= 0
    if np.sum(imag_freq > 0):
        if verbose:
            print(
                f"{np.sum(imag_freq > 0)} imaginary frequencies ignored: {frequencies[imag_freq]}"
            )
        frequencies[imag_freq] = np.nan
    clip_freq = np.clip(frequencies, f_cutoff, np.nan)
    return clip_freq


def calc_zero_point_energy(frequencies):
    """Calculates zero point correction in Hartree/Particle
    Parameters:
    frequencies : array
        Array/List of frequencies in cm^-1
    Returns:
        Zero-Point Correction in Hartree/Particle
    """
    conversion_factor = SPEED_OF_LIGHT * 100
    zpe = np.nansum(0.5 * PLANCK_CONSTANT * frequencies * conversion_factor)
    return zpe / HARTREE2JOULES


def calc_translational_entropy(molar_mass, temperature=298.15, M=None, p=None):
    """Calculates translational component of Entropy
    S_t = R (ln((2pi m kT/h^2)^(3/2) * V) + 5/2)
    with the mass of the molecule m in kg and the volume V in m^3
    following 'Cramer, C. J. Essentials of Computational Chemistry: Theories
    and Models, 2nd ed' p:362, eq:10.18
    Parameters:
    molar_mass : float
        Molar mass of Molecule in atomic units
    temperature : float
        Temperature in Kelvin
    M : float (optional)
        Concentration of Standard State in mol/L
    p : float (optional)
        Pressure of Standard State in atm
    Returns:
        Translational Contribution to Entropy in Cal/Mol-Kelvin
    """

    if p and M:
        raise Warning("Choose Concentration OR Pressure for Standard State")
    elif M:
        V = 1 / (M * 1000 * AVOGARDOS_CONSTANT)  # m^3
    elif p:
        V = GAS_CONSTANT * temperature / (p * ATM2PASCAL * AVOGARDOS_CONSTANT)  # m^3
    else:
        raise ValueError("No Standard State specified")

    return (
        GAS_CONSTANT
        * (
            np.log(
                (
                    (2 * np.pi * molar_mass * AMU2KG * BOLTZMANN_CONSTANT * temperature)
                    / PLANCK_CONSTANT ** 2
                )
                ** (3 / 2)
                * V
            )
            + 2.5
        )
    ) * JOULES2CAL


def calc_vibrational_entropy(frequencies, temperature=298.15):
    """Calculates vibrational component of Entropy
    S_v = R \sum( hv / (kT (exp(hv/kT) - 1)) - ln(1 - exp(-hv/kT)))
    'Cramer, C. J. Essentials of Computational Chemistry: Theories
    and Models, 2nd ed' p:365, eq:10.30
    Parameters:
    frequencies : list
        Frequencies in cm^-1
    temperature : float
        Temperature in Kelvin
    Returns:
        Vibrational Contribution to Entropy in Cal/Mol-Kelvin
    """

    energy_factor = PLANCK_CONSTANT * SPEED_OF_LIGHT * 100
    thermal_energy = BOLTZMANN_CONSTANT * temperature

    energies = frequencies * energy_factor / thermal_energy

    return (
        np.nansum(
            GAS_CONSTANT * energies / (np.exp(energies) - 1)
            - GAS_CONSTANT * np.log(1 - np.exp(-energies))
        )
        * JOULES2CAL
    )


def calc_translational_energy(temperature):
    return 1.5 * GAS_CONSTANT * temperature * JOULES2CAL / 1000


def calc_rotational_energy(temperature):
    return 1.5 * GAS_CONSTANT * temperature * JOULES2CAL / 1000


def calc_vibrational_energy(frequencies, temperature):
    """Calculates Vibrational Energy including ZPE
    U_v = R \sum(hv/2k + hv/k * 1/(exp(hv/kT) - 1))
    Args:
        frequencies (numpy.Array): Frequencies in cm^-1
        temperature (float): Temperature in K
    Returns:
        float: Vibrational Energy in KCal/Mol
    """

    vib_temp = PLANCK_CONSTANT * frequencies * SPEED_OF_LIGHT * 100 / BOLTZMANN_CONSTANT
    vib_energy = GAS_CONSTANT * np.nansum(
        vib_temp / 2 + vib_temp / (np.exp(vib_temp / temperature) - 1)
    )
    return vib_energy * JOULES2CAL / 1000


class Thermochemistry:
    """
    A class to calculate thermodynamical properties from a Gaussian LOG-file
    ...
    Attributes
    ----------
    log_file : str
        Path to LOG-file
    f_cutoff : float (default None)
        Frequency cut off value in cm^-1
    standard_state_M : float (default None)
        Standard state in mol/L
    standard_state_p : float (default None)
        Standard state in atm
    """

    def __init__(
        self,
        log_file,
        f_cutoff=None,
        standard_state_M=None,
        standard_state_p=None,
        temperature=None,
        verbose=True,
    ):
        self.log_file = log_file
        self.f_cutoff = f_cutoff
        self.standard_state_M = standard_state_M
        self.standard_state_p = standard_state_p
        self.verbose = verbose
        self.T = temperature

    def read_lines(self):
        with open(self.log_file, "r") as f:
            lines = f.readlines()
        self.lines = lines

    def normal_termination(self):
        # check if terminated successfully
        if not "Normal termination of Gaussian" in next(
            s for s in reversed(self.lines) if s != "\n"
        ):
            print(f"Abnormal Termination: of {self.log_file}")
            return False
        else:
            return True

    def read_properties(self):
        T0, p0 = get_temperature_and_pressure(self.lines)  # Kelvin
        electronic_energy = get_electronic_energy(self.lines)  # Hartree/Particle
        solv_energy = get_solv_energy(
            self.lines
        )  # kcal/mol THIS IS ONLY THE NON-ELECTROSTATICS!
        frequencies = get_frequencies(self.lines)  # cm^-1
        Sr = get_rotational_entropy(self.lines)  # Cal/Mol-Kelvin
        m = get_molar_mass(self.lines)
        self.T0 = T0
        self.p0 = p0
        self.electronic_energy = electronic_energy
        self.solv_energy = solv_energy
        self.frequencies = frequencies
        self.S_rot = Sr
        self.molar_mass = m
        if not self.T:
            self.T = T0

    def handle_frequencies(self):
        c_freq = clip_frequencies(self.frequencies, self.f_cutoff, verbose=self.verbose)
        self.c_frequencies = c_freq

    def calc_entropies(self):
        St = calc_translational_entropy(
            self.molar_mass, self.T, self.standard_state_M, self.standard_state_p
        )  # Cal/Mol-Kelvin

        if self.verbose:
            stxt = (
                f"{self.standard_state_M} M"
                if self.standard_state_M
                else f"{self.standard_state_p} atm"
            )
            print(f"Calculating the Gibbs Free Energy at {self.T} K")
            print(f"with a Standard State of {stxt}")
            if self.f_cutoff:
                print(f"and with a frequency cutoff of {self.f_cutoff} cm^-1")

        Sv = calc_vibrational_entropy(
            frequencies=self.c_frequencies, temperature=self.T
        )  # Cal/Mol-Kelvin

        S = St + Sv + self.S_rot  # Cal/Mol-Kelvin

        self.S_trans = St
        self.S_vib = Sv
        self.S_tot = S

    def calc_thermal_energies(self):
        trans_energy = calc_translational_energy(self.T)  # KCal/Mol
        rot_energy = calc_rotational_energy(self.T)  # KCal/Mol
        vib_energy = calc_vibrational_energy(self.c_frequencies, self.T)  # KCal/Mol
        tot_energy = np.sum([trans_energy, rot_energy, vib_energy])  # KCal/Mol

        self.U_trans = trans_energy
        self.U_rot = rot_energy
        self.U_vib = vib_energy
        self.U_tot = tot_energy

    def _print(self):
        table = [
            [
                "Electronic Energy ()",
                self.electronic_energy,
                self.electronic_energy * HARTEE2CALMOL / 1000,
            ],
            [
                "Thermal Correction to Gibbs Energy",
                self.gibbs_correction,
                self.gibbs_correction * HARTEE2CALMOL / 1000,
            ],
            [
                "Gibbs Free Energy",
                self.gibbs_free_energy,
                self.gibbs_free_energy * HARTEE2CALMOL / 1000,
            ],
        ]

        print(
            tabulate(
                table,
                headers=["Property", "[Hartree/Particle]", "[KCal/Mol]"],
                numalign="right",
                floatfmt=".4f"
            )
        )

    def run(self):
        """
        Calculates the Gibbs Free Energy in Hartree from a Gaussian LOG-file with the
        option to adjust the Standard State and treat low frequencies as proposed by
        Truhlar and Cramer (doi.org/10.1021/jp205508z, p:14559, bottom right)
        """
        self.read_lines()
        if not self.normal_termination():
            raise Exception("Abnormal Termination of Gaussian")
        self.read_properties()
        self.handle_frequencies()
        self.calc_thermal_energies()
        self.calc_entropies()
        self.zpe = calc_zero_point_energy(self.c_frequencies)

        # Calculate Corrections
        self.thermal_correction_energy = (
            self.U_tot * 1000 / HARTEE2CALMOL
        )  # Hartree/Particle
        self.thermal_correction_enthalpy = (
            self.thermal_correction_energy
            + BOLTZMANN_CONSTANT * self.T / HARTREE2JOULES
        )

        self.entropy_correction = self.S_tot * self.T / HARTEE2CALMOL

        self.gibbs_correction = (
            self.thermal_correction_enthalpy - self.entropy_correction
        )

        # Calculate Enthalpy in Hartree/Particle
        self.enthalpy = self.electronic_energy + self.thermal_correction_enthalpy

        # Calculate Gibbs Free Energy in Hartree/Particle
        self.gibbs_free_energy = (
            self.electronic_energy
            + self.thermal_correction_enthalpy
            - self.entropy_correction
        )

        if self.verbose:
            self._print()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Calculates the Gibbs Free Energy in Hartree."
    )
    parser.add_argument(
        "-f", "--file", help="Gaussian log-file", required=True, type=str
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        help="Frequency cut off value in cm^-1",
        required=False,
        type=float,
    )
    parser.add_argument(
        "-m", "--molarity", help="Standard state in mol/L", required=False, type=float
    )
    parser.add_argument(
        "-p", "--pressure", help="Standard state in atm", required=False, type=float
    )
    args = parser.parse_args()

    thermo = Thermochemistry(
        args.file,
        f_cutoff=args.cutoff,
        standard_state_M=args.molarity,
        standard_state_p=args.pressure,
        verbose=True,
    )
    thermo.run()
