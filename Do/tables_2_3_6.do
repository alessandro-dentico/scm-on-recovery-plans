* ==============================================================================
*   PAPER:  Fiscal Discipline, Uneven Regional Consequences: Evidence from the
*           Italian National Health Service
*           Alessandro Dentico and Cinzia Di Novi
*   FILE:   tables_2_3_6.do
* ==============================================================================
*
*   This file is one part of a larger replication package for the paper
*   above; it is not meant to be run as a standalone script.
*
*   >>> THIS FILE PRODUCES TABLES 2, 3, AND 6 OF THE PAPER <<<
*
*    Table 2  Year-specific estimated treatment effects on the SSN
*             personnel rate
*           	(LaTeX label: tab:effects_personnel)
*
*    Table 3  Year-specific estimated treatment effects on patient
*             out-migration
*             	(LaTeX label: tab:effects_out_migr)
*
*    Table 6  Year-specific placebo-outcome estimates for the total
*             fertility rate
*              	(LaTeX label: tab:effects_fert_rate)
*
*     All three tables share the same row layout: a header row naming
*     every treated region (alphabetical order), one row per
*     post-treatment year, followed by three summary rows per treated
*     region: pre-treatment RMSPE, the joint standardised permutation
*     p-value over the whole post-treatment period, and the same joint
*     p-value restricted to years through 2019 (excluding the COVID
*     period, which complicates causal inference).
*
*     As the code runs, it also prints a results table per region and
*     outcome to the Stata output/log, in red so it stands out from
*     everything else. It is meant to be checked by eye against the
*     tables copied into the paper, so that any accidental change is 
*	  easy to catch.
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
log using "Log/tables_2,3,6", text replace



* ==============================================================================
* 0. SETTING
* ==============================================================================

* ------------------------------------------------------------------------------
* 0.1: define folder paths
* ------------------------------------------------------------------------------

global panel_file "Data/Panel/panel_HFA.dta"  // analysis panel (input)
global table_dir  "Output"                    // root folder for all tables

global treated_order "13 18 15 12 7 14 1 16 20 19"
* Abruzzo Calabria Campania Lazio Liguria Molise Piemonte Puglia Sardegna Sicilia
* Alphabetical order; identical to treated_order in tables_4_5_7.do. Both
* the header row and every data row in write_effects_table step through
* this same list, so the column a region's numbers land in and the
* column its name is printed in can never drift apart.

* ------------------------------------------------------------------------------
* 0.1b: region-code to region-name dictionary (globals, not a value label)
* ------------------------------------------------------------------------------
* Region names are looked up through these globals rather than through
* Stata's own value-label machinery (": label lbl_region <code>"): that
* extended macro function has been observed, from inside a program, to
* silently fall back to the bare numeric code rather than the name, even
* when the label was defined beforehand at file scope. Globals are
* visible inside every program unconditionally and do not depend on
* value-label resolution at all, so they cannot fail in that same way.

global pipername1  "Piemonte"
global pipername2  "Valle d'Aosta"
global pipername3  "Lombardia"
global pipername4  "Trentino-A.A."
global pipername5  "Veneto"
global pipername6  "Friuli-V.G."
global pipername7  "Liguria"
global pipername8  "Emilia-Romagna"
global pipername9  "Toscana"
global pipername10 "Umbria"
global pipername11 "Marche"
global pipername12 "Lazio"
global pipername13 "Abruzzo"
global pipername14 "Molise"
global pipername15 "Campania"
global pipername16 "Puglia"
global pipername17 "Basilicata"
global pipername18 "Calabria"
global pipername19 "Sicilia"
global pipername20 "Sardegna"

* ------------------------------------------------------------------------------
* 0.2: check that the analysis dataset is present
* ------------------------------------------------------------------------------
* "quietly" only suppresses the normal printed output of use/xtdescribe; if
* panel_HFA.dta is missing, Stata still stops with an error here.

quietly {
    use "${panel_file}", clear
    xtdescribe
}
clear



* ==============================================================================
* 1. PROGRAM: run_scm_table
* ==============================================================================
* Builds the synthetic control for one treated region and one outcome,
* twice: once over the whole post-treatment period, and once restricted
* to years through 2019. Posts one row per post-treatment year to the
* open postfile (opened before this program is called, in section 3
* below), and prints a colour-highlighted results table to the Stata
* output/log.
*
* Unlike figures_1-3.do, every command here runs quietly, including
* synth_runner itself: this file runs it twice per region-outcome
* combination purely to extract summary numbers, and the red table
* printed in step 1.5 is the only thing meant to be read.
*
* ARGUMENTS:
*   trunit_num  : numeric region code of the treated region (e.g. 12 for Lazio)
*   outcome_var : name of the outcome variable (e.g. out_migr)
*   first_year  : first year of the analysis window
*   last_year   : last year of the analysis window
*   T_year      : treatment year, i.e. the region's PdR entry year

capture program drop run_scm_table
program define run_scm_table

    args trunit_num outcome_var first_year last_year T_year

    * --------------------------------------------------------------------------
    * 1.1: load the analysis dataset and restrict to the relevant sample
    * --------------------------------------------------------------------------
    * Reloaded fresh here, and again in step 1.4, so that nothing from one
    * step carries over into the other. Wrapped in "quietly" so that this
    * routine bookkeeping does not clutter the output/log; the one line
    * that matters is printed on purpose in step 1.5.

    quietly {
        use "${panel_file}", clear

        keep if year >= `first_year'
        keep if year <= `last_year'
        keep if region_id == `trunit_num' | donor == 1

        xtset region_id year
    }

    local region_label "${pipername`trunit_num'}"
    local T_minus1 = `T_year' - 1

    * --------------------------------------------------------------------------
    * 1.2: build the predictor list passed to synth_runner
    * --------------------------------------------------------------------------
    * Every pre-treatment year of the outcome itself is used as a
    * predictor (see the paper's empirical strategy for why).

    local lag_list ""
    forvalues yr = `first_year'/`T_minus1' {
        local lag_list "`lag_list' `outcome_var'(`yr')"
    }

    * --------------------------------------------------------------------------
    * 1.3: build the synthetic control over the whole post-treatment period
    * --------------------------------------------------------------------------
    * This run covers every post-treatment year through last_year. gen_vars
    * additionally saves each unit's pre-treatment RMSPE as a variable in
    * the dataset (named pre_rmspe), read off below for the treated region.

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

    * --------------------------------------------------------------------------
    * 1.4: build the synthetic control again, restricted to years through 2019
    * --------------------------------------------------------------------------
    * max_lead() tells synth_runner to only consider the first N
    * post-treatment years when computing the joint p-value; N is set here
    * to cover exactly T_year through 2019, so the result excludes the
    * COVID period. Only the resulting joint p-value is needed from this
    * run, so gen_vars is not used.

    quietly {
        use "${panel_file}", clear
        keep if year >= `first_year'
        keep if year <= `last_year'
        keep if region_id == `trunit_num' | donor == 1
        xtset region_id year
    }

    local lead_2019 = 2019 - `T_year' + 1
    assert `lead_2019' >= 1 & `lead_2019' <= `n_leads'

    quietly synth_runner `outcome_var'                                        ///
        `lag_list',                                                          ///
        trunit(`trunit_num') trperiod(`T_year')                              ///
        xperiod(`first_year'(1)`T_minus1')                                   ///
        max_lead(`lead_2019')

    local pval_joint_2019 = e(pval_joint_post_std)

    * --------------------------------------------------------------------------
    * 1.5: print a colour-highlighted results table, and post each row to
    *      the open postfile
    * --------------------------------------------------------------------------
    * "as error" is used purely for its red display style, not because
    * anything went wrong.
    * One display line and one posted row per post-treatment year, so the
    * numbers checked on screen are exactly the numbers written to the
    * table file by write_effects_table. 
	* pre_rmspe, pval_joint_full, and pval_joint_2019 are region-level values 
	* (the same for every year of a given region); they are posted on every row 
    * so that write_effects_table can retrieve them with a simple summarize, 
    * and printed once at the bottom of the table for this region.

    display as error _n "`outcome_var', `region_label'"
    display as error "  Year    Estimate   Std. p-value"

    forvalues l = 1/`n_leads' {
        local yr = `T_year' + `l' - 1

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
* Reads the postfile produced by run_scm_table calls for one outcome, and
* writes effects_table.tex: a header row naming every treated region,
* one row per post-treatment year, then three summary rows (pre-treatment
* RMSPE, joint p-value over the full post-treatment period, joint
* p-value through 2019).
*
* Symbols on the year rows follow the paper's own convention: a dagger
* marks a period-specific standardised p-value of 0 (no donor-region
* placebo produces an equally large effect that year); an asterisk marks
* a value of 0.1 (exactly one placebo does).
*
* ARGUMENTS:
*   wf_path     : path to the postfile .dta to load
*   outcome_var : name of the outcome variable
*   out_dir     : directory where the table file should be saved

capture program drop write_effects_table
program define write_effects_table

    args wf_path outcome_var out_dir

    quietly use "`wf_path'", clear
    quietly levelsof year, local(years)

    * The postfile loaded above is a separate .dta file that never
    * carried over the lbl_region value label from panel_HFA.dta, so it
    * is redefined here, identically to tables_4_5_7.do's
    * write_weight_table, purely so that "label lbl_region <code>" below
    * resolves to a region name (e.g. "Abruzzo") rather than silently
    * falling back to the bare numeric code -- which, left unfixed, also
    * makes every data-row lookup below fail to match region_label and
    * so leaves the whole table blank.

    label define lbl_region         ///
        1  "Piemonte"               ///
        2  "Valle d'Aosta"          ///
        3  "Lombardia"              ///
        4  "Trentino-A.A."          ///
        5  "Veneto"                 ///
        6  "Friuli-V.G."            ///
        7  "Liguria"                ///
        8  "Emilia-Romagna"         ///
        9  "Toscana"                ///
        10 "Umbria"                 ///
        11 "Marche"                 ///
        12 "Lazio"                  ///
        13 "Abruzzo"                ///
        14 "Molise"                 ///
        15 "Campania"               ///
        16 "Puglia"                 ///
        17 "Basilicata"             ///
        18 "Calabria"               ///
        19 "Sicilia"                ///
        20 "Sardegna", replace

    * A literal "$" substituted via a local macro rather than written as
    * "\$" in a file write string: with local macros interpolated in
    * between, "\$" has been observed to leave a stray backslash in the
    * output on the closing delimiter (e.g. "$0.220\$" instead of
    * "$0.220$"), which breaks LaTeX math mode. Substituting a macro that
    * already holds a bare "$" sidesteps that escaping issue entirely.
    local dollar "$"

    * --------------------------------------------------------------------------
    * 2.1: write the LaTeX table
    * --------------------------------------------------------------------------
    * Written as loose rows, not a full table environment, so they can be
    * pasted straight into paper.tex's own table/tabularx wrapper -- but
    * the header row naming every treated region is written here too
    * (not left for paper.tex to supply), so this file is self-contained
    * and the column order can never drift out of sync with whatever is
    * pasted into the paper. Both the header and every data row below
    * step through treated_order (section 0.1), in alphabetical order,
    * rather than any data-derived ordering.

    local out_tex "`out_dir'/effects_table.tex"
    capture file close fh
    file open fh using "`out_tex'", write replace

    file write fh "% Outcome: `outcome_var'" _n

    file write fh "Year"
    foreach t_id of global treated_order {
        local t_label "${pipername`t_id'}"
        file write fh " & \rotatebox{45}{\small `t_label'}"
    }
    file write fh " \\" _n
    file write fh "\midrule" _n

    foreach yr of local years {
        file write fh "`yr'"
        foreach t_id of global treated_order {
            local r "${pipername`t_id'}"
            quietly summarize effect if region_label == "`r'" & year == `yr'
            if r(N) > 0 {
                local est : display %8.3f r(mean)
                local est = trim("`est'")
                quietly summarize pval_std if region_label == "`r'" & year == `yr'
                local pv = r(mean)
                if `pv' == 0          local sym "^{\dagger}"
                else if `pv' <= 0.10  local sym "^{*}"
                else                  local sym ""
                file write fh " & `dollar'`est'`sym'`dollar'"
            }
            else {
                file write fh " &"
            }
        }
        file write fh " \\" _n
    }

    file write fh "\midrule" _n

    file write fh "Pre-treatment RMSPE"
    foreach t_id of global treated_order {
        local r "${pipername`t_id'}"
        quietly summarize pre_rmspe if region_label == "`r'"
        local val : display %8.4f r(mean)
        local val = trim("`val'")
        file write fh " & `val'"
    }
    file write fh " \\" _n

    file write fh "Joint std.\ `dollar'p`dollar'-value (full period)"
    foreach t_id of global treated_order {
        local r "${pipername`t_id'}"
        quietly summarize pval_joint_full if region_label == "`r'"
        local val : display %8.3f r(mean)
        local val = trim("`val'")
        file write fh " & `val'"
    }
    file write fh " \\" _n

    file write fh "Joint std.\ `dollar'p`dollar'-value (through 2019)"
    foreach t_id of global treated_order {
        local r "${pipername`t_id'}"
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
* run_scm_table is called once per treated region, separately for each of
* the three variables, in between opening and closing a postfile for that
* variable. Argument order inside each call:
*   trunit_num  outcome_var  first_year  last_year  T_year
*
* PdR entry years used below: Abruzzo, Campania, Lazio, Molise, Sardegna,
* Sicilia, and Liguria enter in 2007; Calabria enters in 2009; Piemonte
* and Puglia enter in 2010.
* ==============================================================================


* ==============================================================================
* 3A. TABLE 2: SSN PERSONNEL RATE  (personnel, 1994-2022)
* ==============================================================================

local out_dir_pers "${table_dir}/Tables of personnel"
shell mkdir -p "`out_dir_pers'"
tempfile rf_pers

postfile results_handle                                                       ///
    str40 region_label int year double effect double pval_std                 ///
    double pre_rmspe double pval_joint_full double pval_joint_2019            ///
    using `rf_pers', replace

run_scm_table 13 personnel 1994 2022 2007   // Abruzzo
run_scm_table 18 personnel 1994 2022 2009   // Calabria
run_scm_table 15 personnel 1994 2022 2007   // Campania
run_scm_table 12 personnel 1994 2022 2007   // Lazio
run_scm_table  7 personnel 1994 2022 2007   // Liguria
run_scm_table 14 personnel 1994 2022 2007   // Molise
run_scm_table  1 personnel 1994 2022 2010   // Piemonte
run_scm_table 16 personnel 1994 2022 2010   // Puglia
run_scm_table 20 personnel 1994 2022 2007   // Sardegna
run_scm_table 19 personnel 1994 2022 2007   // Sicilia

postclose results_handle
write_effects_table `rf_pers' personnel "`out_dir_pers'"



* ==============================================================================
* 3B. TABLE 3: PATIENT OUT-MIGRATION RATE  (out_migr, 1999-2023)
* ==============================================================================

local out_dir_mig "${table_dir}/Tables of out-migration"
shell mkdir -p "`out_dir_mig'"
tempfile rf_mig

postfile results_handle                                                       ///
    str40 region_label int year double effect double pval_std                 ///
    double pre_rmspe double pval_joint_full double pval_joint_2019            ///
    using `rf_mig', replace

run_scm_table 13 out_migr 1999 2023 2007   // Abruzzo
run_scm_table 18 out_migr 1999 2023 2009   // Calabria
run_scm_table 15 out_migr 1999 2023 2007   // Campania
run_scm_table 12 out_migr 1999 2023 2007   // Lazio
run_scm_table  7 out_migr 1999 2023 2007   // Liguria
run_scm_table 14 out_migr 1999 2023 2007   // Molise
run_scm_table  1 out_migr 1999 2023 2010   // Piemonte
run_scm_table 16 out_migr 1999 2023 2010   // Puglia
run_scm_table 20 out_migr 1999 2023 2007   // Sardegna
run_scm_table 19 out_migr 1999 2023 2007   // Sicilia

postclose results_handle
write_effects_table `rf_mig' out_migr "`out_dir_mig'"



* ==============================================================================
* 3C. TABLE 6: PLACEBO -- TOTAL FERTILITY RATE  (fert_rate, 1980-2024)
* ==============================================================================

local out_dir_fert "${table_dir}/Tables of TFT"
shell mkdir -p "`out_dir_fert'"
tempfile rf_fert

postfile results_handle                                                       ///
    str40 region_label int year double effect double pval_std                 ///
    double pre_rmspe double pval_joint_full double pval_joint_2019            ///
    using `rf_fert', replace

run_scm_table 13 fert_rate 1980 2024 2007   // Abruzzo
run_scm_table 18 fert_rate 1980 2024 2009   // Calabria
run_scm_table 15 fert_rate 1980 2024 2007   // Campania
run_scm_table 12 fert_rate 1980 2024 2007   // Lazio
run_scm_table  7 fert_rate 1980 2024 2007   // Liguria
run_scm_table 14 fert_rate 1980 2024 2007   // Molise
run_scm_table  1 fert_rate 1980 2024 2010   // Piemonte
run_scm_table 16 fert_rate 1980 2024 2010   // Puglia
run_scm_table 20 fert_rate 1980 2024 2007   // Sardegna
run_scm_table 19 fert_rate 1980 2024 2007   // Sicilia

postclose results_handle
write_effects_table `rf_fert' fert_rate "`out_dir_fert'"



* ==============================================================================
* END OF .DO FILE
* ==============================================================================

log close
