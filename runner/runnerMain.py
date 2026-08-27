from runnerClass import *
import re
import numpy as np


def main():
    runner = Runner()
    nml_name1 = "physical_parameters"
    nml_name2 = "calculation_parameters"

    # Bz = [0.01]
    # for b in Bz:
    #     runner.run_slurm_param_value([(nml_name1, "n0_trapped", b)], runsDir="negativeU/nowe/size/80", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 1.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 9.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 8.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 7.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 6.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 5.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 4.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 3.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 2.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    runner.run_slurm_param_value([(nml_name1, "sigma", 10.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 15.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 14.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 13.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 12.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 11.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
if __name__ == "__main__":
    main()