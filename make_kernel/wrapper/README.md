In addition to the compile instructions from Quentin I have edited 
his wrapper to take generic profiles. The profiles are processed in: 
    - ``./rrtmg_lw/run_custom_rrtmg_lw.py``
    - ``./rrtmg_sw/run_custom_rrtmg_sw.py``

These are the only files that need to be edited to process custom 
atmospheric profiles. If other aerosol/chemical profiles are needed to 
be customized, then I will have to edit the wrapper further.

The following files also have a line that needs to be replaced with 
the correct relative directory path:
    - ``./rrtmg_lw/run_custom_rrtmg_lw.py``...... line 8
    - ``./rrtmg_lw/rrtmg_cld_func_band_N.py``.... line 156
    - ``./rrtmg_lw/rrtmg_cld_func_band.py``...... line 106
    - ``./rrtmg_lw/run_rrtmg_cloud_allsky.py``... line 6
    - ``./rrtmg_sw/run_custom_rrtmg_sw.py``...... line 9
    - ``./rrtmg_sw/rrtmg_cld_func_band_N.py``.... line 156
    - ``./rrtmg_sw/rrtmg_clr_func_band.py``...... line 110
    - ``./rrtmg_sw/run_rrtmg_cloud_allsky.py``... line 7

Once those files are customized to use atmospheric profiles, the
script simply needs to be run using ``python run_custom_rrtmg_lw.py``
or ``python run_custom_rrtmg_sw.py`` (when in the correct directories).

For an example of how to create a shell script that automates the 
process of running over several cases, please see below:

(this uses computecanada's sbatch system, see their documentation
for more information)

The text below would be saved as a .sh file and 
run using the command ```sbatch example.sh```

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --account=ctb-itan
#SBATCH --mem=25G
#SBATCH --output=./run_files/write_kernel_lw.out
#SBATCH --error=./run_files/write_kernel_lw.err
#SBATCH --mail-user=catherine.stauffer@mcgill.ca
#SBATCH --mail-type=ALL

module load netcdf;module load python/3.11.5;module load scipy-stack/2023b;source ~/CLS/bin/activate

for m in 01 02 03 04 05 06 07 08 09 10 11 12
do
    python run_custom_rrtmg_lw.py $m
done
