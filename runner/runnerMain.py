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
    runner.run_slurm_param_value([(nml_name1, "sigma", 2.0)], runsDir="STO-self_focusing/tests_focusing", machine="default")
    
if __name__ == "__main__":
    main()