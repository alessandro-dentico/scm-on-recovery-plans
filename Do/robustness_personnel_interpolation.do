* ==============================================================================
*   PAPER:  Fiscal Discipline, Uneven Regional Consequences: Evidence from the
*           Italian National Health Service
*           Alessandro Dentico and Cinzia Di Novi
*   FILE:   robustness_personnel_interpolation.do
* ==============================================================================
*
*   This file is one part of a larger replication package for the paper
*   above; it is not meant to be run as a standalone script.
*
*   >>> THIS FILE IS A ROBUSTNESS CHECK, NOT ONE OF THE PAPER'S NUMBERED
*       TABLES <<<
*
*     data_preparation.do fills a genuine reporting gap in ISTAT indicator
*     9110 (personnel): every region is missing a value for 2014 and 2015,
*     and those two years are filled in by straight-line interpolation
*     between 2013 and 2016 (see data_preparation.do, section 3.2). This
*     file re-runs the personnel SCM specification from tables_2_3_6.do
*     with those two years dropped entirely instead of interpolated, to
*     check whether the interpolation is driving any of the results.
*     Fertility and out-migration have no such gap and are not affected,
*     so this file only covers personnel.
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
log using "Log/robustness_personnel_interpolation", text replace



* ==============================================================================
* 0. SETTING
* ==============================================================================

* ------------------------------------------------------------------------------
* 0.1: define folder paths
* ------------------------------------------------------------------------------

global panel_file "Data/Panel/panel_HFA.dta"  // analysis panel (input)
global table_dir  "Output"                    // root folder for all tables

* ------------------------------------------------------------------------------
* 0.2: check that the analysis dataset is present
* ------------------------------------------------------------------------------

quietly {
    use "${panel_file}", clear
    xtdescribe
}
clear



* ==============================================================================
* 1. PROGRAM: run_scm_robustness
* ==============================================================================
* Builds the synthetic control for one treated region, for personnel
* only, twice: once over the whole post-treatment period, and once
* restricted to years through 2019 -- both with 2014 and 2015 dropped
* from the sample entirely, for every region. Posts one row per
* surviving post-treatment year to the open postfile (opened before this
* program is called, in section 3 below), and prints a colour-highlighted
* results table to the Stata output/log.
*
* Every command here runs quietly, for the same reason as in
* tables_2_3_6.do: the red table printed in step 1.5 is the only thing
* meant to be read.
*
* ARGUMENTS:
*   trunit_num  : numeric region code of the treated region (e.g. 12 for Lazio)
*   outcome_var : name of the outcome variable (always "personnel" here,
*                 but left as an argument for consistency with the other
*                 files in this package)
*   first_year  : first year of the analysis window
*   last_year   : last year of the analysis window
*   T_year      : treatment year, i.e. the region's PdR entry year

capture program drop run_scm_robustness
program define run_scm_robustness

    args trunit_num outcome_var first_year last_year T_year

    * --------------------------------------------------------------------------
    * 1.1: load the analysis dataset, restrict to the relevant sample, and
    *      drop the two interpolated years
    * --------------------------------------------------------------------------
    * drop if inlist(year, 2014, 2015) removes those two years for every
    * region in the sample (treated and donor alike), not just the
    * treated region, so every unit's post-treatment fit is judged on the
    * same set of years.

    quietly {
        use "${panel_file}", clear

        keep if year >= `first_year'
        keep if year <= `last_year'
        keep if region_id == `trunit_num' | donor == 1
        drop if inlist(year, 2014, 2015)

        xtset region_id year
    }

    local region_label : label lbl_region `trunit_num'
    local T_minus1 = `T_year' - 1

    * --------------------------------------------------------------------------
    * 1.2: build the predictor list passed to synth_runner
    * --------------------------------------------------------------------------
    * Every pre-treatment year of the outcome itself is used as a
    * predictor, exactly as in tables_2_3_6.do. 2014 and 2015 are always
    * post-treatment years in this file, so this list is unaffected by
    * dropping them.

    local lag_list ""
    forvalues yr = `first_year'/`T_minus1' {
        local lag_list "`lag_list' `outcome_var'(`yr')"
    }

    * --------------------------------------------------------------------------
    * 1.3: build the synthetic control over the whole (gapped)
    *      post-treatment period
    * --------------------------------------------------------------------------
    * gen_vars saves each unit's pre-treatment RMSPE as a variable in the
    * dataset (named pre_rmspe), read off below for the treated region.

    quietly synth_runner `outcome_var'                                        ///
        `lag_list',                                                          ///
        trunit(`trunit_num') trperiod(`T_year')                              ///
        xperiod(`first_year'(1)`T_minus1')                                   ///
        gen_vars

    matrix B               = e(b)
    matrix PVALS_STD       = e(pvals_std)
    local  n_leads         = colsof(B)
    local  pval_joint_full = e(pval_joint_post_std)

    quietly summarize pre_rmspe if region_id == `trunit_num'
    local pre_rmspe = r(mean)

    * The post-treatment years that actually survived the drop, for the
    * treated region, in calendar order. Column l of B and of PVALS_STD
    * corresponds to the l-th year in this list, not to T_year + l - 1,
    * since 2014 and 2015 are missing from the middle of the sequence.

    quietly levelsof year if region_id == `trunit_num' & year >= `T_year', ///
        local(post_years)
    local n_post : word count `post_years'
    assert `n_leads' == `n_post'

    * --------------------------------------------------------------------------
    * 1.4: build the synthetic control again, restricted to years through
    *      2019
    * --------------------------------------------------------------------------
    * The sample itself is truncated at 2019 here (rather than using
    * synth_runner's max_lead() option, as tables_2_3_6.do does), so that
    * the joint p-value is computed over exactly the surviving years
    * through 2019, regardless of how max_lead() would count periods in a
    * gapped panel. Only the resulting joint p-value is needed from this
    * run, so gen_vars is not used.

    quietly {
        use "${panel_file}", clear
        keep if year >= `first_year'
        keep if year <= 2019
        keep if region_id == `trunit_num' | donor == 1
        drop if inlist(year, 2014, 2015)
        xtset region_id year
    }

    quietly synth_runner `outcome_var'                                        ///
        `lag_list',                                                          ///
        trunit(`trunit_num') trperiod(`T_year')                              ///
        xperiod(`first_year'(1)`T_minus1')

    local pval_joint_2019 = e(pval_joint_post_std)

    * --------------------------------------------------------------------------
    * 1.5: print a colour-highlighted results table, and post each row to
    *      the open postfile
    * --------------------------------------------------------------------------
    * "as error" is used purely for its red display style, not because
    * anything went wrong -- it is the standard way to make this table
    * stand out from the rest of the output in Stata's results window.
    * pre_rmspe, pval_joint_full, and pval_joint_2019 are region-level
    * values (the same for every year of a given region); they are posted
    * on every row so that write_effects_table can retrieve them with a
    * simple summarize, and printed once at the bottom of the table for
    * this region.

    display as error _n "`outcome_var' (no 2014-2015), `region_label'"
    display as error "  Year    Estimate   Std. p-value"

    forvalues l = 1/`n_leads' {
        local yr : word `l' of `post_years'

        display as error %6.0f `yr' "   " %8.3f B[1,`l'] "   " %8.3f PVALS_STD[1,`l']

        post results_handle                                                   ///
            ("`region_label'") (`yr')                                        ///
            (B[1,`l']) (PVALS_STD[1,`l'])                                    ///
            (`pre_rmspe') (`pval_joint_full') (`pval_joint_2019')
    }

    display as error "  Pre-treatment RMSPE:           " %8.4f `pre_rmspe'
    display as error "  Joint p-value (through 2019):   " %6.3f `pval_joint_2019'
    display as error "  Joint p-value (full period):    " %6.3f `pval_joint_full'

end



* ==============================================================================
* 2. PROGRAM: write_effects_table
* ==============================================================================
* Identical in purpose to write_effects_table in tables_2_3_6.do: reads
* the postfile produced by run_scm_robustness calls, and writes
* effects_table.tex -- one row per surviving post-treatment year, then
* three summary rows (pre-treatment RMSPE, joint p-value over the full
* post-treatment period, joint p-value through 2019). There is no row
* for 2014 or 2015, since those years are not in the postfile at all.
*
* ARGUMENTS:
*   wf_path     : path to the postfile .dta to load
*   outcome_var : name of the outcome variable
*   out_dir     : directory where the table file should be saved

capture program drop write_effects_table
program define write_effects_table

    args wf_path outcome_var out_dir

    quietly use "`wf_path'", clear
    quietly levelsof region_label, local(regions)
    quietly levelsof year,         local(years)

    * --------------------------------------------------------------------------
    * 2.1: write the LaTeX table
    * --------------------------------------------------------------------------
    * Written as loose rows, not a full table environment, and not tied
    * to any existing table in paper.tex, since this is a robustness
    * check rather than one of the paper's numbered tables.

    local out_tex "`out_dir'/effects_table.tex"
    capture file close fh
    file open fh using "`out_tex'", write replace

    file write fh "% Outcome: `outcome_var', 2014-2015 excluded (no interpolation)" _n

    foreach yr of local years {
        file write fh "`yr'"
        foreach r of local regions {
            quietly summarize effect if region_label == "`r'" & year == `yr'
            if r(N) > 0 {
                local est : display %8.3f r(mean)
                local est = trim("`est'")
                quietly summarize pval_std if region_label == "`r'" & year == `yr'
                local pv = r(mean)
                if `pv' == 0          local sym "^{\dagger}"
                else if `pv' <= 0.10  local sym "^{*}"
                else                  local sym ""
                file write fh " & \$`est'`sym'\$"
            }
            else {
                file write fh " &"
            }
        }
        file write fh " \\" _n
    }

    file write fh "\midrule" _n

    file write fh "Pre-treatment RMSPE"
    foreach r of local regions {
        quietly summarize pre_rmspe if region_label == "`r'"
        local val : display %8.4f r(mean)
        local val = trim("`val'")
        file write fh " & `val'"
    }
    file write fh " \\" _n

    file write fh "Joint std.\ \$p\$-value (full period)"
    foreach r of local regions {
        quietly summarize pval_joint_full if region_label == "`r'"
        local val : display %8.3f r(mean)
        local val = trim("`val'")
        file write fh " & `val'"
    }
    file write fh " \\" _n

    file write fh "Joint std.\ \$p\$-value (through 2019)"
    foreach r of local regions {
        quietly summarize pval_joint_2019 if region_label == "`r'"
        local val : display %8.3f r(mean)
        local val = trim("`val'")
        file write fh " & `val'"
    }
    file write fh " \\" _n

    file close fh

end



* ==============================================================================
* 3. IMPLEMENTATION
* ==============================================================================
* run_scm_robustness is called once per treated region, for personnel
* only. Argument order inside each call:
*   trunit_num  outcome_var  first_year  last_year  T_year
*
* Same PdR entry years as tables_2_3_6.do: 
* Abruzzo, Campania, Lazio, Molise, Sardegna, Sicilia, and Liguria enter in 2007; 
* Calabria enters in 2009; Piemonte and Puglia enter in 2010.
* ==============================================================================

local out_dir_pers "${table_dir}/Robustness (personnel, no interpolation)"
shell mkdir -p "`out_dir_pers'"
tempfile rf_pers

postfile results_handle                                                       ///
    str40 region_label int year double effect double pval_std                 ///
    double pre_rmspe double pval_joint_full double pval_joint_2019            ///
    using `rf_pers', replace

run_scm_robustness 13 personnel 1994 2022 2007   // Abruzzo
run_scm_robustness 18 personnel 1994 2022 2009   // Calabria
run_scm_robustness 15 personnel 1994 2022 2007   // Campania
run_scm_robustness 12 personnel 1994 2022 2007   // Lazio
run_scm_robustness  7 personnel 1994 2022 2007   // Liguria
run_scm_robustness 14 personnel 1994 2022 2007   // Molise
run_scm_robustness  1 personnel 1994 2022 2010   // Piemonte
run_scm_robustness 16 personnel 1994 2022 2010   // Puglia
run_scm_robustness 20 personnel 1994 2022 2007   // Sardegna
run_scm_robustness 19 personnel 1994 2022 2007   // Sicilia

postclose results_handle
write_effects_table `rf_pers' personnel "`out_dir_pers'"



* ==============================================================================
* END OF .DO FILE
* ==============================================================================

log close
