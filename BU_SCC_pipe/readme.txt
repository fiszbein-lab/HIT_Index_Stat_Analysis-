README
By Xingpei Zhang and Zachary Wakefield
July 16, 2022

For performing differential alternative splicing calculation across phenotype:
Recommended with greater than or equal to 3 replicates for each phenotype

For use with BU SCC:
In HitStat.sh:
1. Chance -P flag, -N flag, depending on project name and selected name for scc run
2. Change cd to the directory with the pipeline
3. Within the HITindex_stat_analysis.py call, change --condition1 and --condition2 to the paths for your selected phenotypes


Otherwise (if running locally or not on BU SCC), follow the following procedure:

cd {path of pipeline}

conda env create -f HITindex_stat_analysis.yml
conda activate HITindex_stat_analysis

python HITindex_stat_analysis.py --statAnalysis --condition1 ./stemcells/cmc1,./stemcells/cmc2,./stemcells/cmc3 --condition2 ./stemcells/psc1,./stemcells/psc2,./stemcells/psc3 --biosignificant 0.3 --minimalExonCount 15 --criticalValue .001 --outlierMethod cooks --outlierDetection 4/n --multipleTesting fdr_bh --output stemcellsDE



More details on options and usage, visit the GitHub: https://github.com/xingpeiz/HITIndex_stat_analysis
