#!/bin/bash

source /etc/profile
module load julia/1.8.5
module load gurobi/gurobi-1000

echo "My SLURM_ARRAY_TASK_ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

julia --project=. fusion_paper/paper_runs/dual_runs/2z_20sc_20periods_dual_100percent_1.10QC_SC/Run.jl $LLSUB_RANK $LLSUB_SIZE 16
