#!/bin/bash
#SBATCH --job-name=mustang_step7
#SBATCH --output=mustang_step7_%j.out
#SBATCH --error=mustang_step7_%j.err
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

source ../../.venv/bin/activate
python deconvolution.py
