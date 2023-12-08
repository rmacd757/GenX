#!/bin/bash

source /etc/profile
module load julia/1.8.5
module load gurobi/gurobi-1000

echo "My SLURM_ARRAY_TASK_ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

julia --project=. fusion_paper/paper_runs/dual_runs/final_report_runs/final_basecase/Run_dual.jl $LLSUB_RANK $LLSUB_SIZE 16
# julia --project=. fusion_paper/paper_runs/dual_runs/final_report_runs/final_uncappedVRE/Run_dual.jl $LLSUB_RANK $LLSUB_SIZE 16
julia --project=. fusion_paper/paper_runs/dual_runs/final_report_runs/final_noNGCCS/Run_dual.jl $LLSUB_RANK $LLSUB_SIZE 16
julia --project=. fusion_paper/paper_runs/dual_runs/final_report_runs/final_baseload/Run_dual.jl $LLSUB_RANK $LLSUB_SIZE 16
julia --project=. fusion_paper/paper_runs/dual_runs/final_report_runs/final_thermstor/Run_dual.jl $LLSUB_RANK $LLSUB_SIZE 16
julia --project=. fusion_paper/paper_runs/dual_runs/final_report_runs/final_lowTrit/Run_dual.jl $LLSUB_RANK $LLSUB_SIZE 16
julia --project=. fusion_paper/paper_runs/dual_runs/final_report_runs/final_withFission/Run_dual.jl $LLSUB_RANK $LLSUB_SIZE 16


