#!/bin/bash
#SBATCH --nodelist=hkbugpusrv[01,02,03]
#SBATCH --nodes=3
#SBATCH --sockets-per-node=2
#SBATCH --cores-per-socket=7
#SBATCH --threads-per-core=2
#SBATCH --ntasks=86
#SBATCH --cpus-per-task=1
#SBATCH --time=60-00:00

echo "Running forecasting_gamma" >> ../progress_record.txt
srun python forecasting_gamma.py

echo "Running forecasting_i5r5" >> ../progress_record.txt
srun python forecasting_i5r5.py

echo "Running forecasting_i5r7" >> ../progress_record.txt
srun python forecasting_i5r7.py

echo "Running forecasting_i7r7" >> ../progress_record.txt
srun python forecasting_i7r7.py

echo "Running forecasting_i6r6" >> ../progress_record.txt
srun python forecasting_i6r6.py

echo "Running forecasting_i7r5" >> ../progress_record.txt
srun python forecasting_i7r5.py

echo "Running forecasting_i5r7_level4" >> ../progress_record.txt
srun python forecasting_i5r7_level4.py

echo "Running forecasting_i5r7_level8" >> ../progress_record.txt
srun python forecasting_i5r7_level8.py

echo "Running forecasting_lognormal" >> ../progress_record.txt
srun python forecasting_lognormal.py

echo "Running forecasting_real" >> ../progress_record.txt
srun python forecasting_real.py

echo "Running vaccination_gamma" >> ../progress_record.txt
srun python vaccination_gamma.py

echo "Running vaccination_i5r5" >> ../progress_record.txt
srun python vaccination_i5r5.py

echo "Running vaccination_i5r7" >> ../progress_record.txt
srun python vaccination_i5r7.py

echo "Running vaccination_i7r7" >> ../progress_record.txt
srun python vaccination_i7r7.py

echo "Running vaccination_lognormal" >> ../progress_record.txt
srun python vaccination_lognormal.py

echo "Running vaccination_real" >> ../progress_record.txt
srun python vaccination_real.py

echo "Running vaccination_i6r6" >> ../progress_record.txt
srun python vaccination_i6r6.py

echo "Running vaccination_i7r5" >> ../progress_record.txt
srun python vaccination_i7r5.py

echo "Running transient_equivalence" >> ../progress_record.txt
srun python transient_equivalence.py

echo "Running calc_r0_g" >> ../progress_record.txt
srun python calc_r0_g.py

echo "Running theory_compare_ijk" >> ../progress_record.txt
srun python theory_compare_ijk.py

echo "Running explain" >> ../progress_record.txt
srun python explain.py

echo "Complete!" >> ../progress_record.txt