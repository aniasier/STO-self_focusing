import f90nml
import textwrap


class RunnerConfig:
    def __init__(self):

        # Using dedent to remove indentation and align every row with the first column of file
        self.jobHeader = {
            "default": textwrap.dedent(
                """\
                #!/bin/bash
                ##### Amount of cores per task
                #SBATCH --cpus-per-task=1
                ##### Partition name
                #SBATCH -p cpu
                ##### Name of job in queuing system
                #SBATCH --job-name=STO
                #SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory
                """
                ),
            "ares": textwrap.dedent(
                """\
                #!/bin/bash -l
                ## Job name
                #SBATCH -J STO
                ## Number of allocated nodes
                #SBATCH -N 1
                ## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
                #SBATCH --ntasks-per-node=1
                ## Memory allocated per core (default is 5GB), comment if mem for whole job should be taken
                ##SBATCH --mem-per-cpu=3800MB
                ## Memory allocated for whole job, comment if mem-per-cpu should be taken
                #SBATCH --mem=16GB
                ## Max task execution time (format is HH:MM:SS)
                #SBATCH --time=72:00:00
                ## Name of grant to which resource usage will be charged
                #SBATCH -A plgqubit-cpu
                ## Name of partition
                #SBATCH -p plgrid
                ## Name of file to which standard output will be redirected
                #SBATCH --output="output.out"
                ## Name of file to which the standard error stream will be redirected
                #SBATCH --error="error.err"
                """
                ),
            "helios": textwrap.dedent(
                """\
                #!/bin/bash -l
                ## Job name
                #SBATCH -J v_tuning_NNN
                ## Number of allocated nodes
                #SBATCH -N 1
                ## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
                #SBATCH --ntasks-per-node=48
                ## Memory allocated per core (default is 5GB), comment if mem for whole job should be taken
                ##SBATCH --mem-per-cpu=3800MB
                ## Memory allocated for whole job, comment if mem-per-cpu should be taken
                #SBATCH --mem=16GB
                ## Max task execution time (format is HH:MM:SS)
                #SBATCH --time=72:00:00
                ## Name of grant to which resource usage will be charged
                #SBATCH -A plgktosto111-cpu
                ## Name of partition
                #SBATCH -p plgrid
                ## Name of file to which standard output will be redirected
                #SBATCH --output="output.out"
                ## Name of file to which the standard error stream will be redirected
                #SBATCH --error="error.err"
                module load GCC/13.2.0 OpenMPI/5.0.3 FlexiBLAS/3.3.1 ScaLAPACK/2.2.0-fb gimkl/2023b
                """
                ),
        }

    def LAO_STO_QD_default_nml(self):
        parser = f90nml.Parser()
        params_nml = parser.reads(
            f"&calculation_parameters \
                Nx=100, \
                Ny=100, \
                Nz=1000, \
                dx=1.0, \
                dz=0.07, \
                MAX_ITER=10000, \
                MAX_ITER_SCF=40, \
                tol=1.0e-8, \
                tol_scf=1.0e-6, \
                alfa=0.5, \
                dt=1e-9, \
                MAX_TIME=500/ \
            &physical_parameters \
                n0_trapped=5.0, \
                L_trapped=10.0, \
                m1=0.28, \
                m2=3.5, \
                norbital=1, \
                sigma=1.0/"
        )
        return params_nml
