#!/bin/bash -l


#$ -P evolution
#$ -l h_rt=99:59:00
#$ -N HITstat_Analysis
#$ -j y
#$ -m ea
#$ -pe omp 28
#$ -l mem_per_core=18G

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

module load miniconda
cd /projectnb2/evolution/zwakefield/HIT_Stat_Analysis

conda env create -f HITindex_stat_analysis.yml
conda activate HITindex_stat_analysis


python HITindex_stat_analysis.py --statAnalysis --condition1 ./stemcells/cmc1,./stemcells/cmc2,./stemcells/cmc3 --condition2 ./stemcells/psc1,./stemcells/psc2,./stemcells/psc3 --biosignificant 0.3 --minimalExonCount 15 --criticalValue .001 --outlierMethod cooks --outlierDetection 4/n --multipleTesting fdr_bh --output stemcellsDE
