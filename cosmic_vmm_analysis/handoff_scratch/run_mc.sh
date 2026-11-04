#!/bin/bash
# usage: run_mc.sh OUT N [extra args]
cd /local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/.claude/worktrees/trigger-acceptance-mc/cosmic_vmm_analysis
export VMM_QA_DIR=/local/home/ak271430/Documents/PostDocSaclay/DAQ_Control_VMM_Beam/vmm_qa
J=/local/home/ak271430/Documents/PostDocSaclay/data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/eff_definitions/json
out=$1; n=$2; shift 2
/usr/bin/time -f "MAXRSS_KB %M  ELAPSED %e s" ~/miniconda3/bin/python trigger_acceptance_mc.py "$out" $J/scintor_run*.json --n "$n" "$@"
