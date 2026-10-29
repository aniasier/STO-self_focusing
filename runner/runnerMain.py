from runnerClass import *
import re
import numpy as np


def main():
    runner = Runner()
    nml_name1 = "physical_parameters"
    nml_name2 = "calculation_parameters"

    Bz = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0]
    for b in Bz:
        runner.run_slurm_param_value([(nml_name1, "sigma", b)], runsDir="STO-self_focusing/oct/sigma_init", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 1.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 9.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 8.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 7.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 6.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 5.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 4.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 3.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 2.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
    # runner.run_slurm_param_value([(nml_name1, "m1", 10.0)], runsDir="STO-self_focusing/epsilon-tests-local/m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 15.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 14.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 13.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 12.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
#     runner.run_slurm_param_value([(nml_name1, "sigma", 11.0)], runsDir="STO-self_focusing/epsilon-changing-newtest/n0-m08", machine="default")
if __name__ == "__main__":
    main()