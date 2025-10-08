#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --account=ctb-itan
#SBATCH --mem=25G
#SBATCH --output=./run_files/write_kernel_lw.out
#SBATCH --error=./run_files/write_kernel_lw.err
#SBATCH --mail-user=catherine.stauffer@mcgill.ca
#SBATCH --mail-type=ALL

module load netcdf;module load python/3.11.5;module load scipy-stack/2023b;source ~/CLS/bin/activate

for m in 1 2 3 4 5 6 7 8 9 10 11 12
do
python run_custom_rrtmg_lw.py $m 160 538
python run_custom_rrtmg_lw.py $m 1450 1499
python run_custom_rrtmg_lw.py $m 2740 509
python run_custom_rrtmg_lw.py $m 1450 1612
python run_custom_rrtmg_lw.py $m 2740 3479
python run_custom_rrtmg_lw.py $m 4030 741
done

# python run_custom_rrtmg_lw.py $m 160 538
# python run_custom_rrtmg_lw.py $m 1670 1285
# python run_custom_rrtmg_lw.py $m 3180 353
# python run_custom_rrtmg_lw.py $m 1670 1612
# python run_custom_rrtmg_lw.py $m 3180 3120
# python run_custom_rrtmg_lw.py $m 4690 593
