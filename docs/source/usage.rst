=====
Usage
=====

To use locuaz in a project::

    import locuaz



SLURM
^^^^^

An example script with SLURM:

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
    source /m100/home/userexternal/pbarlett/.bashrc
    conda activate locuaz
    module load profile/lifesc
    module load autoload gromacs/2021.4

    python /home/user/locuaz/locuaz/locuaz.py config.yaml


PBS
^^^^

An example script with PBS:

.. code-block:: console

    #!/bin/bash
    #PBS -N locuaz
    #PBS -l walltime=00:15:00
    #PBS -l select=1:ncpus=20:ngpus=2:mpiprocs=20
    #PBS -q debug

    cd $PBS_O_WORKDIR
    export OMP_NUM_THREADS=4
    source /home/pbarletta/.bashrc
    module load gromacs/2021.4
    module load mpi
    conda activate locuaz

    python /home/user/locuaz/locuaz/locuaz.py config.yaml
