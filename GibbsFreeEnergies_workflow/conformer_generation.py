#!/groups/kemi/mharris/.conda/envs/rdkit_2020_09/bin/python

import sys
import os
import subprocess
import re

import pandas as pd
import numpy as np

from string import ascii_lowercase as alc


from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem


def embed_mult_conf(smiles, nconfs=None):
    ps = AllChem.ETKDGv3()
    mol = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(mol)
    if not nconfs:
        nrot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        nconfs = 3+3*nrot
        print(nconfs)

    mol = Chem.AddHs(mol)
    rdmolops.AssignStereochemistry(mol)
    confs = AllChem.EmbedMultipleConfs(mol, numConfs=nconfs, params=ps)
    return mol, charge


def write_confs_xyz(mol, molID):
    file_names = []
    for confid in range(len(mol.GetConformers())):
        file_name = str(molID)+'_'+str(confid)+'.xyz'
        file_names.append(file_name)
        rdmolfiles.MolToXYZFile(mol, file_name, confId=confid)
    return file_names


def run_cmd(cmd):
    """
    Run command line
    """
    cmd = cmd.split()
    print(cmd)
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output, err = p.communicate()
    return output.decode('utf-8')



def do_xtb_calculations(file_names, charge, mol_id):
    energies = [] #hartree
    lowest_e_conformer = None
    lowest_e = None
    with open(str(mol_id)+"_conformers.xyz", "w") as mol_file:
        for xyz_file in file_names:
            output = run_cmd("/groups/kemi/koerstz/opt/xtb/6.1/bin/xtb {0} --opt tight --gbsa methanol --gfn2 --chrg {1}".format(xyz_file, charge))
            with open(xyz_file[:-4]+'_gfn2.log', 'w') as _file:
                _file.write(output)
            with open("xtbopt.xyz",'r') as fi:
                #mol_file.write(fi.read())
                line = fi.readline()
                while line:
                    if "SCF done" in line:
                        e = np.float(line.split()[2])
                        print(e)
                    line = fi.readline()
            energies.append(e)
            if (not lowest_e) or (e < lowest_e):
                lowest_e = e
                lowest_e_conformer = xyz_file
            os.rename('xtbopt.xyz', xyz_file[:-4]+'_opt.xyz')
    os.rename(lowest_e_conformer, str(mol_id)+"_lowest_xtb_conf.xyz")
    return energies, str(mol_id)+"_lowest_xtb_conf.xyz"


def write_sp_com_file(xyz_file, cpus, mem, charge):
    """ prepares com file for gaussian """

    com_file = xyz_file[:-4]+'_opt.com'
    with open(com_file, 'w') as _file:
        _file.write('%nprocshared='+str(cpus)+'\n')
        _file.write('%mem='+str(mem)+'\n')
        _file.write('#B3LYP/6-31+G(d,p) scrf=(smd, solvent=methanol) empiricaldispersion=gd3 \n\n')
        _file.write('something title\n\n')
        _file.write('{0} 1\n'.format(charge))
        with open(xyz_file, 'r') as file_in:
            lines = file_in.readlines()[2:]
            _file.writelines(lines)
        _file.write('\n')

    return com_file




def write_opt_com_file(xyz_file, cpus, mem, charge):
    """ prepares com file for gaussian """

    com_file = xyz_file[:-4]+'_opt.com'
    with open(com_file, 'w') as _file:
        _file.write('%nprocshared='+str(cpus)+'\n')
        _file.write('%mem='+str(mem)+'\n')
        _file.write('#opt freq B3LYP/6-31+G(d,p) scrf=(smd, solvent=methanol) empiricaldispersion=gd3 \n\n')
        _file.write('something title\n\n')
        _file.write('{0} 1\n'.format(charge))
        with open(xyz_file, 'r') as file_in:
            lines = file_in.readlines()[2:]
            _file.writelines(lines)
        _file.write('\n')

    return com_file


def calc_gaussian(com_file_writer, *args):
    """ Do gaussian-xtb ts optimization """
    com_file = com_file_writer(*args)
    output = run_cmd("srun /software/kemi/g16C/g16 {0}".format(com_file))
    with open(com_file[:-4]+'.out', 'w') as _file:
        _file.write(output)

    return com_file[:-4]+'.out'





if __name__ == "__main__":
    os.environ["XTBHOME"] = "/groups/kemi/koerstz/opt/xtb/6.1/bin"
    os.environ["OMP_STACKSIZE"] = '8G'
    os.environ["OMP_NUM_THREADS"] = '4'
    os.environ["MKL_NUM_THREADS"] = '4'
    os.system('ulimit -s unlimited')

    #df = pd.read_csv(sys.argv[1], index_col=0)

    mol_id = sys.argv[1]
    SMILES = sys.argv[2]
    cpus = sys.argv[3]
    mem = sys.argv[4]

    os.mkdir(str(mol_id))
    os.chdir(str(mol_id))
    mol, charge = embed_mult_conf(SMILES)
    xyz_files = write_confs_xyz(mol, mol_id)
    energies, xyz_file = do_xtb_calculations(xyz_files, charge, mol_id)
    out_file = calc_gaussian(write_opt_com_file, xyz_file, cpus, mem, charge)
    os.chdir("../")






