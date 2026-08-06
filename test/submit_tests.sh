#!/bin/bash
#SBATCH --job-name=nf-test-suite
#SBATCH --partition=hsph
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=/n/home04/smaharjan/biobakery-nextflow/test/results/slurm_%j.out
#SBATCH --error=/n/home04/smaharjan/biobakery-nextflow/test/results/slurm_%j.err

set -euo pipefail

# Load environment
source /n/lab_storage/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
module use /n/lab_storage/huttenhower_lab/tools/hutlab/src/modules_rocky8
module load jdk/21.0.2-fasrc01

# Confirm nextflow is reachable
/n/lab_storage/huttenhower_lab/tools/nextflow/24.10.4/bin/nextflow -version

# Run the test suite from the repo root
cd /n/home04/smaharjan/biobakery-nextflow
bash test/run_tests.sh
