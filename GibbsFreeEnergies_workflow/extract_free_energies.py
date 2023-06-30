import tarfile
import sys
import os
import shutil
import thermochemistry
import xyz2mol_local
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import PeriodicTable
from rdkit.Chem import GetPeriodicTable

def get_free_energy(dft_file, standard_state=1, freq_cutoff=20):

    thermo = thermochemistry.Thermochemistry(dft_file, f_cutoff=freq_cutoff,
                                             standard_state_p=standard_state,
                                             verbose=False)

    thermo.run()

    return thermo.gibbs_free_energy


def get_sp_energy(dft_file):
    """
    Extracts the single point electronic energy from Gaussian ouput
    """
    scf_energy = None
    with open(dft_file, 'r') as _file:
        line = _file.readline()
        while line:
            if 'SCF Done' in line:
                scf_energy = np.float(line.split()[4])
            line = _file.readline()

    return scf_energy


def extract_optimized_structure(out_file, n_atoms, atom_labels):
    """
    After waiting for the constrained optimization to finish, the
    resulting structure from the constrained optimization is
    extracted and saved as .xyz file ready for TS optimization.
    """
    optimized_xyz_file = out_file[:-4]+".xyz"
    with open(out_file, 'r') as ofile:
        line = ofile.readline()
        while line:
            if 'Standard orientation' in line or 'Input orientation' in line:
                coordinates = np.zeros((n_atoms, 3))
                for i in range(5):
                    line = ofile.readline()
                for i in range(n_atoms):
                    coordinates[i, :] = np.array(line.split()[-3:])
                    line = ofile.readline()
            line = ofile.readline()
    with open(optimized_xyz_file, 'w') as _file:
        _file.write(str(n_atoms)+'\n\n')
        for i in range(n_atoms):
            _file.write(atom_labels[i])
            for j in range(3):
                _file.write(' '+"{:.5f}".format(coordinates[i, j]))
            _file.write('\n')

    return optimized_xyz_file


def check_acs(xyz_file_in, gaus_file_out, charge, use_huckel=True):

    pt = GetPeriodicTable()
    atoms, _, xyz_coordinates = xyz2mol_local.read_xyz_file(xyz_file_in)
    atom_labels = []
    for a in atoms:
        atom_labels.append(PeriodicTable.GetElementSymbol(pt, a))
    xyz_file_out = extract_optimized_structure(gaus_file_out, len(atoms),
                                               atom_labels)
    atoms2, _, xyz_coordinates2 = xyz2mol_local.read_xyz_file(xyz_file_out)
    try:
        AC1, _ = xyz2mol_local.xyz2AC(atoms, xyz_coordinates, charge,
                                      use_huckel=use_huckel)
        AC2, _ = xyz2mol_local.xyz2AC(atoms2, xyz_coordinates2, charge,
                                      use_huckel=use_huckel)
    except:
        return "Error: could not get AC from xyz file"
    if np.all(AC1 == AC2):
        return None

    return "Error: AC changes during DFT optimization"


if __name__ == "__main__":

    pwd = os.getcwd()
    input_csv = sys.argv[1]
    df = pd.read_csv(input_csv, index_col=0)
    for idx, smiles in zip(df.index, df.can_smiles):
        mol = Chem.MolFromSmiles(smiles)
        charge = Chem.GetFormalCharge(mol)
        os.mkdir(str(idx))
        try:
            tar = tarfile.open(str(idx)+'.tar.gz', 'r:gz')
            tar.extractall()
            tar.close()

        except:
            df.loc[idx, "error_message"] = "Error: no space left on device"
            shutil.rmtree(str(idx))
            continue
        os.chdir(str(idx))
        dft_file = str(idx)+"_lowest_xtb_conf_opt.out"
        xyz_file_in = str(idx)+"_lowest_xtb_conf.xyz"
        try:
            ac_check = check_acs(xyz_file_in, dft_file, charge)
            print(ac_check)
            df.loc[idx, "error_message"] = str(ac_check)
            e_Ha = get_free_energy(dft_file)
            print(idx, e_Ha)
            df.loc[idx, 'Gibbs_Ha'] = e_Ha
            e_sp_hartree = get_sp_energy(dft_file)
            print(idx, e_sp_hartree)
            df.loc[idx, 'sp_hartree'] = e_sp_hartree
        except:
            if not os.path.exists(xyz_file_in):
                df.loc[idx, "error_message"] = "Error: can't embed"
            else:
                df.loc[idx, "error_message"] = "Error: something wrong in Gaussian opt"
        os.chdir("../")
        shutil.rmtree(str(idx))
    df.to_csv("test_df.csv")
    print(df)
