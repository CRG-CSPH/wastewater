#!/bin/bash
#SBATCH --output=ww_mobility_download.log
#SBATCH --error=ww_mobility_download.err
#SBATCH --exclude=csphbiostats.ucdenver.pvt

source /home/hillandr/venv/bin/activate

PYTHONUNBUFFERED=1 python3 download.py
