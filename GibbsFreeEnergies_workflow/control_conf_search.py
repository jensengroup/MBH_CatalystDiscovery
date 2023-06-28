#!/groups/kemi/mharris/.conda/envs/rdkit_2020_09/bin/python


import os
import textwrap
import sys
import time
import random

import pandas as pd



def qsub_prep(script_path, cpus, mem, smiles_idx, smiles):
    """
    write qsub file for SLURM subsmissin
    """
    pwd = os.getcwd()

    qsub_file = """\
    #!/bin/sh
    #SBATCH --job-name={3}
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task={1}
    #SBATCH --mem={2}
    #SBATCH --error={5}/{3}.stderr
    #SBATCH --output={5}/{3}.stdout
    #SBATCH --ntasks=1
    #SBATCH --time=1000:00:00
    #SBATCH --partition=kemi1
    #SBATCH --no-requeue

    #mkdir /scratch/$SLURM_JOB_ID
    #cp ../initial_structures.xyz /scratch/$SLURM_JOB_ID

    export GAUSS_SCRDIR=/scratch/$SLURM_JOB_ID

    cd /scratch/$SLURM_JOB_ID

    #run python code

    export GAUSS_SCRDIR=/scratch/$SLURM_JOB_ID

    ({0} {3} '{4}' {1} {2})
    #cp output back

    tar -zcvf {3}.tar.gz {3}

    cp {3}.tar.gz {5}/{3}.tar.gz

    rm {5}/{3}_qsub.tmp

    #rm -r /scratch/$SLURM_JOB_ID

    """.format(script_path, cpus, mem, smiles_idx, smiles, pwd)

    with open(str(smiles_idx) + "_qsub.tmp", "w") as qsub:
        qsub.write(textwrap.dedent(qsub_file))

    return str(smiles_idx) + "_qsub.tmp"


def run_calculations(smiles_df, script_path, max_queue, cpus, mem):
    """
    For each given smiles - submit a conformational search. Only submit new jobs when less than max_queue jobs in the queue
    """
    submitted_jobs = set()
    for idx in smiles_df.index:
        smiles = smiles_df.loc[idx, 'can_smiles']
        #if idx not in submitted_smiles:
        #    os.mkdir(str(idx))
        #    submitted_smiles.add(idx)
        qsub_name = qsub_prep(script_path, cpus, mem, idx, smiles)
        slurmid = os.popen("sbatch " + qsub_name).read()
        slurmid = int(slurmid.strip().split()[-1])

        submitted_jobs.add(slurmid)

        if len(submitted_jobs) >= max_queue:
            while True:
                job_info = os.popen("squeue -u mharris").readlines()[1:]
                current_jobs = {int(job.split()[0]) for job in job_info}

                if len(current_jobs) >= max_queue:
                    time.sleep(30)
                else:
                    finished_jobs = submitted_jobs - current_jobs
                    print("finished jobs: ", finished_jobs)
                    for job in finished_jobs:
                        submitted_jobs.remove(job)
                    break


if __name__ == "__main__":
    SMILES_LIST = sys.argv[1]
    SMILES_DF = pd.read_csv(SMILES_LIST, index_col=0)

    CPUS = 4
    MEM = "8GB"
    MAX_QUEUE = 195

    SCRIPT = '/groups/kemi/mharris/github/get_free_energies/conformer_generation.py'

    _DIR = "conformer_calcs"

    os.mkdir(_DIR)
    os.chdir(_DIR)

    run_calculations(SMILES_DF, SCRIPT, MAX_QUEUE, CPUS, MEM)

