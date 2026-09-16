* ==============================================================================
*   PAPER:  Fiscal Discipline, Uneven Regional Consequences: Evidence from the
*           Italian National Health Service
*           Alessandro Dentico and Cinzia Di Novi
*   FILE:   figures_1-3.do
* ==============================================================================
*
*   This file is one part of a larger replication package for the paper
*   above; it is not meant to be run as a standalone script.
*
*   >>> THIS FILE PRODUCES FIGURES 1-3 OF THE PAPER <<<
*
*    Figure 1  Observed and synthetic public healthcare personnel rates
*           	(LaTeX label: fig:personnel_trajectories)
*
*    Figure 2  Observed and synthetic patient out-migration rates
*              	(LaTeX label: fig:outmigration_trajectories)
*
*    Figure 3  Placebo synthetic controls for the total fertility rate
*              	(LaTeX label: fig:tft_placebo)
*
*     Figure 1 shows the PdR's effect on the SSN personnel rate. It
*     contains 10 panels, (a)-(j), one per treated region.
*
*     Figure 2 shows the PdR's effect on the patient out-migration rate.
*     It contains 10 panels, (a)-(j), one per treated region.
*
*     Figure 3 shows the placebo on the total fertility rate. It
*     contains 10 panels, (a)-(j), one per treated region.
*
*     For each outcome and each treated region, the code below produces
*     exactly one such panel: a PDF plotting that region's observed
*     outcome over time against its synthetic twin, built through the
*     Synthetic Control Method (SCM) as described in the paper. Each panel
*     is saved as Output/Graphs of <outcome>/<region>.pdf. For example,
*     the Abruzzo panel of Figure 1 is Output/Graphs of personnel/abruzzo.pdf.
*
* ==============================================================================



* ------------------------------------------------------------------------------
* Preliminary commands
* ------------------------------------------------------------------------------

clear all
set more off
version 17

* $drive is set in master.do, which calls this file; it points to the
* replication package's root folder.
cd "$drive"

capture log close
log using "Log/figures_1-3", text replace



* ==============================================================================
* 0. SETTING
* ==============================================================================

* ------------------------------------------------------------------------------
* 0.1: define folder paths
* ------------------------------------------------------------------------------

global panel_file "Data/Panel/panel_HFA.dta"  // analysis panel (input)
global graph_dir  "Output"                    // root folder for all graphs

* ------------------------------------------------------------------------------
* 0.2: check that the analysis dataset is present
* ------------------------------------------------------------------------------
* panel_HFA.dta is produced by data_preparation.do. This file only reads it,
* so this is  just a quick check that it exists and has the expected panel shape
* before running anything else.

use "${panel_file}", clear
xtdescribe
clear



* ==============================================================================
* 1. PROGRAM: run_scm_graph
* ==============================================================================
* Builds one synthetic control and saves one treated-vs-synthetic graph.
* Called once per treated region, separately for each of the three
* variables (personnel, out_migr, fert_rate) in section 2 below.
*
* ARGUMENTS:
*   trunit_num  : numeric region code of the treated region (e.g. 12 for Lazio)
*   outcome_var : name of the outcome variable (e.g. out_migr)
*   first_year  : first year of the analysis window
*   last_year   : last year of the analysis window
*   T_year      : treatment year, i.e. the region's PdR entry year
*   T_end       : last year of the region's PdR; pass "" if the region is
*                 still under the PdR at the end of the analysis window
*   graph_subdir: output subfolder name under Output/ for this variable
*                 (e.g. "Graphs of personnel")
*   ytitle_str  : y-axis label to print on the graph
*   y_min       : lower bound for the y-axis; pass "" for an automatic scale
*   y_max       : upper bound for the y-axis; pass "" for an automatic scale
*                 Supply both to lock the same scale across all regions for
*                 one variable, which makes panels visually comparable.
*
* OUTPUT:
*   One PDF per call, saved as Output/<graph_subdir>/<region name>.pdf.

capture program drop run_scm_graph
program define run_scm_graph

    args trunit_num outcome_var first_year last_year T_year T_end graph_subdir ytitle_str y_min y_max

    * --------------------------------------------------------------------------
    * 1.1: load the analysis dataset and restrict to the relevant sample
    * --------------------------------------------------------------------------
    * Reloaded fresh on every call so that changes made for one region (e.g.
    * keeping only certain years) never carry over to the next call.

    use "${panel_file}", clear

    keep if year >= `first_year'
    keep if year <= `last_year'

    * Keep only the treated region itself and the never-treated donor pool.
    * Every other treated region is dropped, so it cannot be used (even
    * accidentally) to build this region's synthetic control.
    keep if region_id == `trunit_num' | donor == 1

    xtset region_id year

    * --------------------------------------------------------------------------
    * 1.2: region name, and where the graph will be saved
    * --------------------------------------------------------------------------
    * region_label reads the region's name from the value label attached in
    * data_preparation.do (e.g. "Abruzzo"). region_file is the same name in
    * lower case, used only for the saved file name (e.g. "abruzzo").

    local region_label : label lbl_region `trunit_num'
    local region_file   = lower("`region_label'")

    local out_path "${graph_dir}/`graph_subdir'"
    shell mkdir -p "`out_path'"

    local T_minus1 = `T_year' - 1

    * --------------------------------------------------------------------------
    * 1.3: build the predictor list passed to synth_runner
    * --------------------------------------------------------------------------
    * This specification uses every pre-treatment year of the outcome
    * itself as a predictor, so the resulting local looks like, e.g.:
    * "out_migr(1999) out_migr(2000) ... out_migr(2006)"

    local lag_list ""
    forvalues yr = `first_year'/`T_minus1' {
        local lag_list "`lag_list' `outcome_var'(`yr')"
    }

    * --------------------------------------------------------------------------
    * 1.4: build the synthetic control
    * --------------------------------------------------------------------------
    * synth_runner builds the region's synthetic control and, via gen_vars, 
	* saves the resulting counterfactual trajectory as a new variable, 
	* `outcome_var'_synth, which section 1.5 plots. 

    synth_runner `outcome_var'                                               ///
        `lag_list',                                                          ///
        trunit(`trunit_num') trperiod(`T_year')                              ///
        xperiod(`first_year'(1)`T_minus1')                                   ///
        gen_vars

    * --------------------------------------------------------------------------
    * 1.5: draw the treated-vs-synthetic graph
    * --------------------------------------------------------------------------
    * Conventions used for every one of the 30 panels, so they all look
    * consistent side by side in the paper: greyscale only (solid black =
    * treated, dashed black = synthetic); legend reads "Treated" /
    * "Synthetic"; x-axis ticks at the first year, last year, PdR entry
    * and exit years, and every multiple of 5 in between; a vertical line
    * marks PdR entry, and a second one marks PdR exit for regions that
    * left before the end of the window (Liguria, Sardegna, Piemonte).

    * -- y-axis title --
    local ylab "`ytitle_str'"

    * -- solid line at PdR exit (if the region left the PdR within the window) --
    local xline_end ""
    if "`T_end'" != "" {
        local xline_end "xline(`T_end', lcolor(gs8) lpattern(solid) lwidth(thin))"
    }

    * -- x-axis tick labels: first year, last year, PdR entry year, PdR exit
    *    year (if any), plus every multiple of 5 strictly in between --
    local xlab_full "`first_year' `T_year' `last_year'"
    if "`T_end'" != "" local xlab_full "`xlab_full' `T_end'"
    forvalues yr = 1980(5)2030 {
        if `yr' > `first_year' & `yr' < `last_year' {
            local xlab_full "`xlab_full' `yr'"
        }
    }

    * -- y-axis range: fixed if both bounds are supplied, automatic
    *    otherwise --
    if "`y_min'" != "" & "`y_max'" != "" {
        local yrange "yscale(range(`y_min' `y_max')) ylabel(`y_min'(`= (`y_max'-`y_min')/4')`y_max', angle(0))"
    }
    else {
        local yrange "ylabel(, angle(0))"
    }

    * -- draw the graph --
    twoway                                                                        ///
        (line `outcome_var'       year if region_id == `trunit_num',             ///
            lcolor(black) lpattern(solid) lwidth(medthick))                      ///
        (line `outcome_var'_synth year if region_id == `trunit_num',             ///
            lcolor(black) lpattern(dash)  lwidth(medthick)),                     ///
        xline(`T_year', lcolor(gs8) lpattern(solid) lwidth(thin))               ///
        `xline_end'                                                               ///
        xlabel(`xlab_full', angle(45))                                           ///
        `yrange'                                                                  ///
        xtitle("Year")                                                            ///
        ytitle("`ylab'", size(medsmall))                                          ///
        legend(order(1 "Treated" 2 "Synthetic")                                  ///
               position(6) rows(1)                                                ///
               region(lcolor(black)) bmargin(tiny) size(medsmall))               ///
        graphregion(color(white)) plotregion(lcolor(none))                        ///
        scheme(s2mono)

    graph export "`out_path'/`region_file'.pdf", replace

end



* ==============================================================================
* 2. IMPLEMENTATION
* ==============================================================================
* run_scm_graph is called once per treated region, separately for each of
* the three variables. Because the layout is identical across variables,
* only the arguments differ from one block to the next.
*
* Argument order inside each call:
*   trunit_num  outcome_var  first_year  last_year  T_year  T_end
*   graph_subdir  ytitle_str  y_min  y_max
*
* PdR entry and exit years used below (see data_preparation.do, section
* 5.1, for how these were established):
*   Abruzzo, Campania, Lazio, Molise, Sicilia   : enter 2007, still under
*                                                  the PdR at end of window
*   Liguria, Sardegna                           : enter 2007, exit 2009
*                                                  (T_end = 2010)
*   Calabria                                    : enter 2009, still under
*                                                  the PdR at end of window
*   Piemonte                                    : enter 2010, exit 2016
*                                                  (T_end = 2017)
*   Puglia                                      : enter 2010, still under
*                                                  the PdR at end of window
* ==============================================================================


* ==============================================================================
* 2A. OUTCOME 1: SSN PERSONNEL RATE  (personnel, 1994-2022)
* ==============================================================================
* Produces Figure 1 (10 panels total: fig:personnel_trajectories). Y-axis
* locked to [0, 20] for every region, so the panels are directly
* comparable to one another.

* 2A.1: Abruzzo (region_id = 13) -- PdR from 2007, ongoing throughout the window
run_scm_graph 13 personnel 1994 2022 2007 "" "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.2: Calabria (region_id = 18) -- PdR from 2009, ongoing throughout the window
run_scm_graph 18 personnel 1994 2022 2009 "" "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.3: Campania (region_id = 15) -- PdR from 2007, ongoing throughout the window
run_scm_graph 15 personnel 1994 2022 2007 "" "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.4: Lazio (region_id = 12) -- PdR from 2007, ongoing throughout the window
run_scm_graph 12 personnel 1994 2022 2007 "" "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.5: Liguria (region_id = 7) -- PdR 2007-2009, exits 2010
run_scm_graph 7 personnel 1994 2022 2007 2010 "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.6: Molise (region_id = 14) -- PdR from 2007, ongoing throughout the window
run_scm_graph 14 personnel 1994 2022 2007 "" "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.7: Piemonte (region_id = 1) -- PdR 2010-2016, exits 2017
run_scm_graph 1 personnel 1994 2022 2010 2017 "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.8: Puglia (region_id = 16) -- PdR from 2010, ongoing throughout the window
run_scm_graph 16 personnel 1994 2022 2010 "" "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.9: Sardegna (region_id = 20) -- PdR 2007-2009, exits 2010
run_scm_graph 20 personnel 1994 2022 2007 2010 "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20

* 2A.10: Sicilia (region_id = 19) -- PdR from 2007, ongoing throughout the window
run_scm_graph 19 personnel 1994 2022 2007 "" "Graphs of personnel" "SSN personnel (per 1,000 inhabitants)" 0 20



* ==============================================================================
* 2B. OUTCOME 2: PATIENT OUT-MIGRATION RATE  (out_migr, 1999-2023)
* ==============================================================================
* Produces Figure 2 (10 panels total: fig:outmigration_trajectories).
* Y-axis bounds vary by region: a single common scale would make the
* trajectories for regions with a low out-migration rate almost flat and
* hard to read, so each region uses the same range of 20 units but with
* different bounds.

* 2B.1: Abruzzo (region_id = 13) -- PdR from 2007, ongoing throughout the window
run_scm_graph 13 out_migr 1999 2023 2007 "" "Graphs of out-migration" "Patient out-migration (%)" 0 20

* 2B.2: Calabria (region_id = 18) -- PdR from 2009, ongoing throughout the window
run_scm_graph 18 out_migr 1999 2023 2009 "" "Graphs of out-migration" "Patient out-migration (%)" 10 30

* 2B.3: Campania (region_id = 15) -- PdR from 2007, ongoing throughout the window
run_scm_graph 15 out_migr 1999 2023 2007 "" "Graphs of out-migration" "Patient out-migration (%)" 0 20

* 2B.4: Lazio (region_id = 12) -- PdR from 2007, ongoing throughout the window
run_scm_graph 12 out_migr 1999 2023 2007 "" "Graphs of out-migration" "Patient out-migration (%)" 0 20

* 2B.5: Liguria (region_id = 7) -- PdR 2007-2009, exits 2010
run_scm_graph 7 out_migr 1999 2023 2007 2010 "Graphs of out-migration" "Patient out-migration (%)" 0 20

* 2B.6: Molise (region_id = 14) -- PdR from 2007, ongoing throughout the window
run_scm_graph 14 out_migr 1999 2023 2007 "" "Graphs of out-migration" "Patient out-migration (%)" 15 35

* 2B.7: Piemonte (region_id = 1) -- PdR 2010-2016, exits 2017
run_scm_graph 1 out_migr 1999 2023 2010 2017 "Graphs of out-migration" "Patient out-migration (%)" 0 20

* 2B.8: Puglia (region_id = 16) -- PdR from 2010, ongoing throughout the window
run_scm_graph 16 out_migr 1999 2023 2010 "" "Graphs of out-migration" "Patient out-migration (%)" 0 20

* 2B.9: Sardegna (region_id = 20) -- PdR 2007-2009, exits 2010
run_scm_graph 20 out_migr 1999 2023 2007 2010 "Graphs of out-migration" "Patient out-migration (%)" 0 20

* 2B.10: Sicilia (region_id = 19) -- PdR from 2007, ongoing throughout the window
run_scm_graph 19 out_migr 1999 2023 2007 "" "Graphs of out-migration" "Patient out-migration (%)" 0 20



* ==============================================================================
* 2C. PLACEBO: TOTAL FERTILITY RATE  (fert_rate, 1980-2024)
* ==============================================================================
* Produces Figure 3 (10 panels total: fig:tft_placebo). The same entry
* and exit years as sections 2A-2B are used.
* Y-axis locked to [0.5, 2.5] for every region.

* 2C.1: Abruzzo (region_id = 13) -- PdR from 2007, ongoing throughout the window
run_scm_graph 13 fert_rate 1980 2024 2007 "" "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.2: Calabria (region_id = 18) -- PdR from 2009, ongoing throughout the window
run_scm_graph 18 fert_rate 1980 2024 2009 "" "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.3: Campania (region_id = 15) -- PdR from 2007, ongoing throughout the window
run_scm_graph 15 fert_rate 1980 2024 2007 "" "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.4: Lazio (region_id = 12) -- PdR from 2007, ongoing throughout the window
run_scm_graph 12 fert_rate 1980 2024 2007 "" "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.5: Liguria (region_id = 7) -- PdR 2007-2009, exits 2010
run_scm_graph 7 fert_rate 1980 2024 2007 2010 "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.6: Molise (region_id = 14) -- PdR from 2007, ongoing throughout the window
run_scm_graph 14 fert_rate 1980 2024 2007 "" "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.7: Piemonte (region_id = 1) -- PdR 2010-2016, exits 2017
run_scm_graph 1 fert_rate 1980 2024 2010 2017 "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.8: Puglia (region_id = 16) -- PdR from 2010, ongoing throughout the window
run_scm_graph 16 fert_rate 1980 2024 2010 "" "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.9: Sardegna (region_id = 20) -- PdR 2007-2009, exits 2010
run_scm_graph 20 fert_rate 1980 2024 2007 2010 "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5

* 2C.10: Sicilia (region_id = 19) -- PdR from 2007, ongoing throughout the window
run_scm_graph 19 fert_rate 1980 2024 2007 "" "Graphs of TFT" "Total fertility rate (children per woman)" 0.5 2.5



* ==============================================================================
* END OF .DO FILE
* ==============================================================================

log close
