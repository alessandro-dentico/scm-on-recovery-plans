* ==============================================================================
*   PAPER:  Fiscal Discipline, Uneven Regional Consequences: Evidence from the
*           Italian National Health Service
*           Alessandro Dentico and Cinzia Di Novi
*   FILE:   tables_4_5_7.do
* ==============================================================================
*
*   This file is one part of a larger replication package for the paper
*   above; it is not meant to be run as a standalone script.
*
*   >>> THIS FILE PRODUCES TABLES 4, 5, AND 7 OF THE PAPER <<<
*
*    Table 4  Synthetic-control donor weights for the personnel rate
*           	(LaTeX label: tab:weights_personnel)
*
*    Table 5  Synthetic-control donor weights for the patient
*             out-migration rate
*             	(LaTeX label: tab:weights_out_migr)
*
*    Table 7  Synthetic-control donor weights for the total fertility
*             rate placebo
*              	(LaTeX label: tab:weights_fert_rate)
*
*     All three tables share the same layout: a header row naming every
*     donor region, one row per treated region (also named), one column
*     per donor region, cell = that donor's weight in the region's
*     synthetic control (blank if zero). Both dimensions run in
*     alphabetical order: Abruzzo, Calabria, Campania, Lazio, Liguria,
*     Molise, Piemonte, Puglia, Sardegna, Sicilia (rows); Basilicata,
*     Emilia-Romagna, Friuli-V.G., Lombardia, Marche, Toscana,
*     Trentino-A.A., Umbria, Valle d'Aosta, Veneto (columns).
*
*     As the code runs, it also prints a results table per region and
*     outcome to the Stata output/log, in red so it stands out from
*     everything else. It is meant to be checked by eye against the
*     tables copied into the paper, so that any accidental change is
*     easy to catch.
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
log using "Log/tables_4,5,7", text replace



* ==============================================================================
* 0. SETTING
* ==============================================================================

* ------------------------------------------------------------------------------
* 0.1: define folder paths and the row/column order used in paper.tex
* ------------------------------------------------------------------------------

global panel_file "Data/Panel/panel_HFA.dta"  // analysis panel (input)
global table_dir  "Output"                    // root folder for all tables

global treated_order "13 18 15 12 7 14 1 16 20 19"
* Abruzzo Calabria Campania Lazio Liguria Molise Piemonte Puglia Sardegna Sicilia

global donor_order "17 8 6 3 11 9 4 10 2 5"
* Basilicata Emilia-Romagna Friuli-V.G. Lombardia Marche Toscana
* Trentino-A.A. Umbria "Valle d'Aosta" Veneto

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

quietly {
    use "${panel_file}", clear
    xtdescribe
}
clear



* ==============================================================================
* 1. PROGRAM: run_scm_weights
* ==============================================================================
* Builds the synthetic control for one treated region and one outcome,
* using the bare synth command (not synth_runner, since only the donor
* weights are needed here, not the permutation-based effect estimates).
* Posts one row per donor to the open postfile (opened before this
* program is called, in section 3 below), and prints a colour-highlighted
* breakdown of the region's synthetic twin to the Stata output/log.
*
* Every command here runs quietly, for the same reason as in this 
* file's sibling, tables_2_3_6.do: the red table printed
* in step 1.5 is the only thing meant to be read.
*
* Donor weights depend only on the pre-treatment window (first_year to
* T_year-1); last_year is passed only so the sample restriction matches
* the rest of the replication package and has no effect on the result.
*
* ARGUMENTS:
*   trunit_num  : numeric region code of the treated region (e.g. 12 for Lazio)
*   outcome_var : name of the outcome variable (e.g. out_migr)
*   first_year  : first year of the analysis window
*   last_year   : last year of the analysis window
*   T_year      : treatment year, i.e. the region's PdR entry year

capture program drop run_scm_weights
program define run_scm_weights

    args trunit_num outcome_var first_year last_year T_year

    * --------------------------------------------------------------------------
    * 1.1: load the analysis dataset and restrict to the relevant sample
    * --------------------------------------------------------------------------

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
    * 1.2: build the predictor list passed to synth
    * --------------------------------------------------------------------------
    * Every pre-treatment year of the outcome itself is used as a
    * predictor, identical to the specification in tables_2_3_6.do.

    local lag_list ""
    forvalues yr = `first_year'/`T_minus1' {
        local lag_list "`lag_list' `outcome_var'(`yr')"
    }

    * --------------------------------------------------------------------------
    * 1.3: build the synthetic control and read off the donor weights
    * --------------------------------------------------------------------------
    * e(W_weights) is a matrix with one row per donor: column 1 is the
    * donor's numeric region code, column 2 is its weight.

    quietly synth `outcome_var' `lag_list',                                    ///
        trunit(`trunit_num') trperiod(`T_year')                               ///
        xperiod(`first_year'(1)`T_minus1')

    matrix W       = e(W_weights)
    local  n_donors = rowsof(W)

    * --------------------------------------------------------------------------
    * 1.4: post each donor's weight to the open postfile, and record it
    *      under its own region code for step 1.5
    * --------------------------------------------------------------------------
    * Posted in whatever order the W_weights matrix lists donors in.
	* wt_<code> is a local macro per donor code, used in step 1.5
    * to print the donors back out in the fixed order defined in
    * section 0.1, regardless of the order synth returned them in.

    local wsum = 0
    forvalues r = 1/`n_donors' {
        local d_id  = W[`r', 1]
        local d_wt  = W[`r', 2]
        local wt_`d_id' = `d_wt'
        local wsum  = `wsum' + `d_wt'

        local d_label "${pipername`d_id'}"
        post results_handle ("`region_label'") ("`d_label'") (`d_wt')
    }

    * --------------------------------------------------------------------------
    * 1.5: print a colour-highlighted breakdown, in the same donor order
    *      used throughout the rest of this file
    * --------------------------------------------------------------------------
    * "as error" is used purely for its red display style, not because
    * anything went wrong. Only donors with a non-negligible weight
    * (> 0.0005) are printed, since those are what actually compose the
    * region's synthetic twin.

    display as error _n "`outcome_var', `region_label'"
    display as error "  Donor               Weight"

    foreach d_id of global donor_order {
        local d_wt = `wt_`d_id''
        if `d_wt' > 0.0005 {
            local d_label "${pipername`d_id'}"
            display as error %-18s "`d_label'" %8.3f `d_wt'
        }
    }

    display as error "  (weights sum to" %7.3f `wsum' ")"

end



* ==============================================================================
* 2. PROGRAM: write_weight_table
* ==============================================================================
* Reads the postfile produced by run_scm_weights calls for one outcome,
* and writes weights_table.tex: one row per treated region, one column
* per donor region.
*
* ARGUMENTS:
*   wf_path     : path to the postfile .dta to load
*   outcome_var : name of the outcome variable
*   out_dir     : directory where the table file should be saved

capture program drop write_weight_table
program define write_weight_table

    args wf_path outcome_var out_dir

    quietly use "`wf_path'", clear

    * --------------------------------------------------------------------------
    * 2.1: write the LaTeX table
    * --------------------------------------------------------------------------
    * The header row and every data row are written by stepping through
    * treated_order and donor_order (section 0.1, both alphabetical),
    * converting each region code to its name via the pipername globals
    * (section 0.1b) as it is used, rather than sorting the region names
    * found in the data.

    local out_tex "`out_dir'/weights_table.tex"
    capture file close fh
    file open fh using "`out_tex'", write replace

    file write fh "% Outcome: `outcome_var'" _n

    * Header row naming every donor region, so the table is self-contained
    * and does not rely on paper.tex separately supplying (and keeping in
    * sync) the same column order. "Region" labels the row-header column;
    * donor columns then follow donor_order (section 0.1), alphabetical,
    * the same list every data row below steps through.

    file write fh "Region"
    foreach d_id of global donor_order {
        local d_label "${pipername`d_id'}"
        file write fh " & `d_label'"
    }
    file write fh " \\" _n
    file write fh "\midrule" _n

    foreach t_id of global treated_order {
        local t_label "${pipername`t_id'}"
        file write fh "`t_label'"
        foreach d_id of global donor_order {
            local d_label "${pipername`d_id'}"
            quietly summarize weight if region_label == "`t_label'" & donor_label == "`d_label'"
            if r(mean) > 0.0005 {
                local val : display %6.3f r(mean)
                local val = trim("`val'")
                file write fh " & \textbf{`val'}"
            }
            else {
                file write fh " &"
            }
        }
        file write fh " \\" _n
    }

    file close fh

end



* ==============================================================================
* 3. IMPLEMENTATION
* ==============================================================================
* run_scm_weights is called once per treated region, separately for each
* of the three variables, in between opening and closing a postfile for
* that variable. Argument order inside each call:
*   trunit_num  outcome_var  first_year  last_year  T_year
*
* PdR entry years used below: Abruzzo, Campania, Lazio, Molise, Sardegna,
* Sicilia, and Liguria enter in 2007; Calabria enters in 2009; Piemonte
* and Puglia enter in 2010.
* ==============================================================================


* ==============================================================================
* 3A. TABLE 4: SSN PERSONNEL RATE  (personnel, 1994-2022)
* ==============================================================================

local out_dir_pers "${table_dir}/Tables of personnel"
shell mkdir -p "`out_dir_pers'"
tempfile wf_pers

postfile results_handle                                                       ///
    str40 region_label str40 donor_label double weight                        ///
    using `wf_pers', replace

run_scm_weights 13 personnel 1994 2022 2007   // Abruzzo
run_scm_weights 18 personnel 1994 2022 2009   // Calabria
run_scm_weights 15 personnel 1994 2022 2007   // Campania
run_scm_weights 12 personnel 1994 2022 2007   // Lazio
run_scm_weights  7 personnel 1994 2022 2007   // Liguria
run_scm_weights 14 personnel 1994 2022 2007   // Molise
run_scm_weights  1 personnel 1994 2022 2010   // Piemonte
run_scm_weights 16 personnel 1994 2022 2010   // Puglia
run_scm_weights 20 personnel 1994 2022 2007   // Sardegna
run_scm_weights 19 personnel 1994 2022 2007   // Sicilia

postclose results_handle
write_weight_table `wf_pers' personnel "`out_dir_pers'"



* ==============================================================================
* 3B. TABLE 5: PATIENT OUT-MIGRATION RATE  (out_migr, 1999-2023)
* ==============================================================================

local out_dir_mig "${table_dir}/Tables of out-migration"
shell mkdir -p "`out_dir_mig'"
tempfile wf_mig

postfile results_handle                                                       ///
    str40 region_label str40 donor_label double weight                        ///
    using `wf_mig', replace

run_scm_weights 13 out_migr 1999 2023 2007   // Abruzzo
run_scm_weights 18 out_migr 1999 2023 2009   // Calabria
run_scm_weights 15 out_migr 1999 2023 2007   // Campania
run_scm_weights 12 out_migr 1999 2023 2007   // Lazio
run_scm_weights  7 out_migr 1999 2023 2007   // Liguria
run_scm_weights 14 out_migr 1999 2023 2007   // Molise
run_scm_weights  1 out_migr 1999 2023 2010   // Piemonte
run_scm_weights 16 out_migr 1999 2023 2010   // Puglia
run_scm_weights 20 out_migr 1999 2023 2007   // Sardegna
run_scm_weights 19 out_migr 1999 2023 2007   // Sicilia

postclose results_handle
write_weight_table `wf_mig' out_migr "`out_dir_mig'"



* ==============================================================================
* 3C. TABLE 7: PLACEBO -- TOTAL FERTILITY RATE  (fert_rate, 1980-2024)
* ==============================================================================

local out_dir_fert "${table_dir}/Tables of TFT"
shell mkdir -p "`out_dir_fert'"
tempfile wf_fert

postfile results_handle                                                       ///
    str40 region_label str40 donor_label double weight                        ///
    using `wf_fert', replace

run_scm_weights 13 fert_rate 1980 2024 2007   // Abruzzo
run_scm_weights 18 fert_rate 1980 2024 2009   // Calabria
run_scm_weights 15 fert_rate 1980 2024 2007   // Campania
run_scm_weights 12 fert_rate 1980 2024 2007   // Lazio
run_scm_weights  7 fert_rate 1980 2024 2007   // Liguria
run_scm_weights 14 fert_rate 1980 2024 2007   // Molise
run_scm_weights  1 fert_rate 1980 2024 2010   // Piemonte
run_scm_weights 16 fert_rate 1980 2024 2010   // Puglia
run_scm_weights 20 fert_rate 1980 2024 2007   // Sardegna
run_scm_weights 19 fert_rate 1980 2024 2007   // Sicilia

postclose results_handle
write_weight_table `wf_fert' fert_rate "`out_dir_fert'"



* ==============================================================================
* END OF .DO FILE
* ==============================================================================

log close
