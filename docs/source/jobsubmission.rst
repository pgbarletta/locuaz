===============
Job submission
===============

While locuaz may be ran on a PC, it was developed to run on UNIX-based clusters where multiple GPUs are available,
and since these usually also include a workload manager, here are 2 sample submission scripts for the 2 most popular
workload managers.


Running from SLURM
^^^^^^^^^^^^^^^^^^^^

Here's an example script with SLURM:

.. code-block:: console

    #!/bin/bash
    #SBATCH -N1
    #SBATCH -n4
    #SBATCH --cpus-per-task=32
    #SBATCH --gres=gpu:4
    #SBATCH --time=24:00:00
    #SBATCH --job-name locuaz
    #SBATCH -o salida_locuaz
    #SBATCH -e error_locuaz
    #SBATCH --exclusive

    cd $SLURM_SUBMIT_DIR
    module load profile/lifesc
    module load autoload gromacs/2021.4
    source /m100/home/userexternal/pbarlett/.bashrc
    conda activate locuaz

    locuaz config.yaml


Running from PBS
^^^^^^^^^^^^^^^^^^

And another one with PBS:

.. code-block:: console

    #!/bin/bash
    #PBS -N locuaz
    #PBS -l walltime=00:15:00
    #PBS -l select=1:ncpus=20:ngpus=2:mpiprocs=20
    #PBS -q debug

    cd $PBS_O_WORKDIR
    export OMP_NUM_THREADS=4
    module load gromacs/2021.4
    module load mpi
    source /home/pbarletta/.bashrc
    conda activate locuaz

    locuaz.py config.yaml
