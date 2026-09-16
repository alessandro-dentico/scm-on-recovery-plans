* ==============================================================================
*   PAPER:  Fiscal Discipline, Uneven Regional Consequences: Evidence from the
*           Italian National Health Service
*           Alessandro Dentico and Cinzia Di Novi
*   FILE:   master.do
* ==============================================================================
*
*   This is the master file for the replication package: running it start
*   to finish reproduces every table and figure in the paper.
*
*   Written and run with Stata 17 on macOS.
*
* ==============================================================================



* ==============================================================================
* BASIC SETUP
* ==============================================================================

clear all
set more off
cap log close
version 17

* ------------------------------------------------------------------------------
* ROOT FOLDER STRUCTURE
* ------------------------------------------------------------------------------
* $drive must point to a folder containing this master.do file plus four
* subfolders:
*   Do/       every other .do file (data_preparation.do, figures_1-3.do, ...)
*   Data/     raw and processed data
*   Output/   every figure and table produced by the pipeline
*   Log/      one text log per .do file

global drive ""   // <-- type the replication package's root folder path here
cd "$drive"

* ------------------------------------------------------------------------------
* Install required packages
* ------------------------------------------------------------------------------
* "cap" (capture) before each install means Stata does not stop with an
* error if a package is already installed.

cap ssc install synth,        all
cap ssc install synth_runner, all



* ==============================================================================
* DATA PREPARATION
* ==============================================================================
* Builds panel_HFA.dta, the analysis panel used by every file below.

do "$drive/Do/data_preparation.do"



* ==============================================================================
* FIGURES 1-3: SCM TRAJECTORIES (PERSONNEL, OUT-MIGRATION, TFR PLACEBO)
* ==============================================================================
* Figure 1: SSN personnel rate
* Figure 2: patient out-migration rate
* Figure 3: total fertility rate placebo

do "$drive/Do/figures_1-3.do"



* ==============================================================================
* TABLES 2, 3, 6: SCM EFFECTS (PERSONNEL, OUT-MIGRATION, TFR PLACEBO)
* ==============================================================================
* Table 2: year-specific effects, SSN personnel rate
* Table 3: year-specific effects, patient out-migration rate
* Table 6: year-specific placebo estimates, total fertility rate

do "$drive/Do/tables_2_3_6.do"



* ==============================================================================
* TABLES 4, 5, 7: SCM DONOR WEIGHTS (PERSONNEL, OUT-MIGRATION, TFR PLACEBO)
* ==============================================================================
* Table 4: donor weights, SSN personnel rate
* Table 5: donor weights, patient out-migration rate
* Table 7: donor weights, total fertility rate placebo

do "$drive/Do/tables_4_5_7.do"



* ==============================================================================
* ROBUSTNESS CHECK: PERSONNEL WITHOUT 2014-2015 INTERPOLATION
* ==============================================================================
* Repeats the personnel SCM specification from tables_2_3_6.do with 2014
* and 2015 dropped instead of interpolated, to check whether the
* interpolation drives any of the results. Not one of the paper's
* numbered tables.

do "$drive/Do/robustness_personnel_interpolation.do"



* ==============================================================================
* END OF MASTER FILE
* ==============================================================================
