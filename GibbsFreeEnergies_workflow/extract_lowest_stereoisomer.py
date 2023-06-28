import sys
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdmolops




def remove_stereo_chemistry_energies(df):
    for i in df.index:
        #print(i)
        mp = Chem.MolFromSmiles(df.loc[i, "can_smiles"])
        rdmolops.RemoveStereochemistry(mp)
        sp = Chem.MolToSmiles(mp)
        df.loc[i, 'smiles_nostereo'] = Chem.MolToSmiles(Chem.MolFromSmiles(sp))

    return df


if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], index_col=0)
    df = remove_stereo_chemistry_energies(df)
    print(df)
    df_min = df.sort_values("sp_hartree").groupby("smiles_nostereo", as_index=False).first()
    #df_min = df.groupby('smiles_nostereo')["sp_hartree", "can_smiles"].min()
    #print(len(df_min))
    #print(len(df.groupby('smiles_nostereo')["sp_hartree"].min()))
    print(df_min)
    df_min.to_csv("minimum_stereo_isomer.csv")
    df_min[["can_smiles"]].to_csv("important_fragment_isomers.csv")
