#!/bin/bash
#SBATCH --job-name=pseudospots
#SBATCH --output=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/pseudospots_generation.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:10:00

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/pseudospots_generation.py