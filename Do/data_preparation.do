* ==============================================================================
*   PAPER:  Fiscal Discipline, Uneven Regional Consequences: Evidence from the
*           Italian National Health Service
*           Alessandro Dentico and Cinzia Di Novi
*   FILE:   data_preparation.do
* ==============================================================================
*
*   This file is one part of a larger replication package for the paper
*   above; it is not meant to be run as a standalone script.
*
*   >>> THIS FILE PRODUCES panel_HFA.dta, THE ANALYSIS PANEL USED BY EVERY
*       OTHER FILE IN THIS REPLICATION PACKAGE <<<
*
*    personnel   SSN personnel rate, per 1,000 inhabitants (ISTAT 9110)
*           	[OUTCOME 1]
*
*    out_migr    Patient out-migration rate, % (ISTAT 7343)
*              	[OUTCOME 2]
*
*    fert_rate   Total fertility rate, children per woman (ISTAT 0180)
*              	[PLACEBO]
*
*    pdr, donor  Treatment (=1 under an active PdR) and donor-pool (=1 if
*              	never treated) indicators, built by hand in this file
*
*     For each of the three raw variables, the code cleans the raw ISTAT
*     file and reshapes it to long (region x year) format. The three
*     cleaned variables are then merged into a single panel, a known gap
*     in personnel is filled by interpolation, two unit conversions are
*     applied, and the treatment and donor-pool indicators are built by
*     hand. The result is saved as Data/Panel/panel_HFA.dta.
*
* ------------------------------------------------------------------------------
* WHERE THE RAW DATA COMES FROM
* ------------------------------------------------------------------------------
* SOURCE
*   All three variables come from ISTAT's "Health for All" (HFA) database,
*   the Italian National Institute of Statistics' repository of regional
*   health indicators.
*
* STEP 1
*   Download the HFA application (only for Windows) from
*   https://www.istat.it/wp-content/uploads/2026/04/HFA_giugno.zip
*
* STEP 2
*   Query each indicator separately in the HFA application, selecting all
*   20 Italian regions and the full available year range; every HFA
*   indicator has a 4-digit code.
*
* STEP 3
*   Export as a .txt table, with one row per region and one column per
*   year.
*
* STEP 4
*   Re-save the exported .txt file as a .csv file, encoded as UTF-8 and
*   using a semicolon (;) as the field separator (a comma would clash
*   with ISTAT's own use of a comma as the decimal separator),
*	and name it after its 4-digit ISTAT code (e.g. 9110.csv). 
*	The row-by-row layout of these files is documented in section 1, 
*	where the import program reads them.
*
* ------------------------------------------------------------------------------
* DATA AVAILABILITY
* ------------------------------------------------------------------------------
*    personnel   1994-2022  annual; gap at 2014-2015 -> interpolated
*    out_migr    1999-2023  annual, no gaps
*    fert_rate   1980-2024  annual, no gaps
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
log using "Log/data_preparation", text replace



* ==============================================================================
* 0. SETTING
* ==============================================================================

* ------------------------------------------------------------------------------
* 0.1: define folder paths
* ------------------------------------------------------------------------------
* Three subfolders inside Data/, one per global below:
* 	- Raw HFA	   : the three raw CSV files described in the header, named
* 	  				 by ISTAT 4-digit indicator code (9110.csv, 7343.csv,
* 	  				 0180.csv)
* 	- Reshaped HFA : intermediate .dta files produced by this file, one per
* 	  				 variable (e.g. out_migr.dta); safe to delete and
*					 regenerate by re-running this file
* 	- Panel	   	   : the final merged panel (panel_HFA.dta)
* Each global below is just a shortcut for its folder path, so that the
* rest of the file can write, e.g., "${raw_dir}/9110.csv" instead of
* typing out the full path every time.

global raw_dir      "Data/Raw HFA"       // raw CSV files
global reshaped_dir "Data/Reshaped HFA"  // per-variable reshaped .dta files
global out_dir      "Data/Panel"         // final master panel



* ==============================================================================
* 1. IMPORT PROGRAM
* ==============================================================================
* Every raw HFA .csv file (regardless of which variable it holds) has an
* identical structure, since every indicator is exported through the same
* process (see the header for how the .csv files were obtained):
*   Row 1     : the indicator's title, written out in full          -> SKIP
*   Row 2     : column headers 										-> SKIP
*   Rows 3-22 : one row per Italian region, in ISTAT's fixed 
*               numeric order (1 = Piemonte, 2 = Valle d'Aosta,
*               ..., 20 = Sardegna), with one column per year
*   Rows 23+  : national and macro-area aggregates (e.g. "Italia",
*               "Nord-ovest") rather than individual regions        -> DROP
* Missing or not-yet-available cells are always left blank, never coded
* as text (e.g. "n.d.") or with a numeric missing-value code.
*
* It is because this layout never changes from one indicator to
* the next that the cleaning steps below are written ONCE, as a small
* program, rather than repeated separately for each variable. Section 2
* simply re-runs this same program on a different raw file for each of
* the three variables.
*
* Arguments:
*   varname  : meaningful column name to assign to the cleaned data
*              (e.g. out_migr)
*   filename : ISTAT code matching the raw CSV filename (e.g. 7343, so the
*              program reads Data/Raw HFA/7343.csv)

capture program drop import_wide_csv
program define import_wide_csv

    args varname filename

    * --------------------------------------------------------------------------
    * 1.1: import the raw CSV as text
    * --------------------------------------------------------------------------
    * varnames(nonames) : the file has no usable header row, so Stata just
    *                      calls the columns v1, v2, v3, ...
    * stringcols(_all)  : import every column as text, not as a number. This
    *                      matters because the cells still contain commas as
    *                      decimal points and text like "n.d." for missing
    *                      values, neither of which Stata could read as a
    *                      number yet.

    import delimited "${raw_dir}/`filename'.csv", ///
        delimiter(";")                             ///
        varnames(nonames)                          ///
        stringcols(_all)                           ///
        encoding("UTF-8")                          ///
        clear

    * --------------------------------------------------------------------------
    * 1.2: drop the title row and the column-header row
    * --------------------------------------------------------------------------
    * Row 1 = indicator title, Row 2 = column headers. Neither contains data.

    drop if _n <= 2

    * --------------------------------------------------------------------------
    * 1.3: keep only the 20 Italian regions (codes 1-20)
    * --------------------------------------------------------------------------
    * Column v1 holds the region code as text at this point. Trimming
    * removes stray spaces, and destring turns it into an actual number so
    * it can be compared to 20. Rows above code 20 are the
    * national/macro-area aggregate rows described above, so they are
    * dropped.

    quietly replace v1 = strtrim(v1)
    destring v1, replace force
    drop if missing(v1) | v1 > 20

    * --------------------------------------------------------------------------
    * 1.4: rename v1 to region_id and attach region-name labels
    * --------------------------------------------------------------------------
    * A value label just tells Stata to display "Piemonte" instead of "1"
    * when it prints region_id, without changing the underlying numeric
    * code.

    rename v1 region_id

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

    label values region_id lbl_region
    label var   region_id "ISTAT region code (1-20)"

    * --------------------------------------------------------------------------
    * 1.5: drop unneeded columns
    * --------------------------------------------------------------------------
    * v2 : region name written out as text, redundant now that region_id
    *      is labelled.
    * v3 : an empty column that always reads "2" or blank in the raw file.
    *
    * NOTE: years 1980-1989 (columns v4-v13) are kept even though personnel
    * and out_migr only start later, so that all three variables go through
    * an identical column-cleaning step below. Years with no data simply
    * end up missing, which does not affect any later analysis.

    drop v2 v3

    * --------------------------------------------------------------------------
    * 1.6: clean and rename the year columns (1980-2025)
    * --------------------------------------------------------------------------
    * The raw file has one column per year, running from 1980 (column v4)
    * to 2025 (column v49); column number = year - 1976. For each year
    * column this loop:
    *   - trims stray spaces
    *   - replaces the Italian decimal comma with a normal decimal point
    *   - turns any blank cell into a true Stata missing value
    *   - destrings: converts the cleaned text into a number, and turns
    *     anything that is still not a valid number (like "n.d.") into
    *     missing
    *   - renames the column from v<col> to y<year>, e.g. v18 becomes y1994

    forvalues yr = 1980/2025 {
        local col = `yr' - 1976
        quietly replace v`col' = strtrim(v`col')
        quietly replace v`col' = subinstr(v`col', ",", ".", 1)
        quietly replace v`col' = "" if regexm(v`col', "^[ ]*$")
        destring v`col', replace force
        rename  v`col' y`yr'
    }

    * --------------------------------------------------------------------------
    * 1.7: reshape from wide to long format
    * --------------------------------------------------------------------------
    * Right now the data is "wide": one row per region, with a separate
    * column for every year (y1980, y1981, ...). "reshape long" turns this
    * into "long" format: one row per region-year pair, with a single
    * column called y holding the value for that region and year. This
    * produces 20 regions x 46 years = 920 rows. The line after renames y
    * to the variable's actual name (e.g. out_migr), passed in as varname.

    reshape long y, i(region_id) j(year)
    rename  y `varname'
    label var `varname' "`varname'"

    * --------------------------------------------------------------------------
    * 1.8: sort and save as an intermediate file
    * --------------------------------------------------------------------------
    * Saved on its own so it can be merged with the other variables in
    * section 3, and so it can be inspected or regenerated independently.

    sort region_id year
    save "${reshaped_dir}/`varname'.dta", replace

end



* ==============================================================================
* 2. IMPORT ALL VARIABLES
* ==============================================================================
* Syntax:  import_wide_csv  varname  filename
*
* Each line below is one call to the program defined in section 1. Because
* every raw file shares the same layout, each call runs through exactly the
* same eight steps (1.1-1.8), and only the input file and the resulting
* variable name differ.

* Outcome 1: SSN personnel rate. Confirmed available 1994-2022 annually,
* with a gap at 2014-2015 for every region (fixed in section 3.2).
import_wide_csv personnel 9110

* Outcome 2: patient out-migration rate. Confirmed available 1999-2023 annually, 
* with no gaps.
import_wide_csv out_migr  7343

* Placebo (fake-outcome) variable: total fertility rate. Not an outcome of
* interest: imported solely for the fake-outcome placebo test. Confirmed
* available 1980-2024 annually, with no known systematic gaps.
import_wide_csv fert_rate 0180



* ==============================================================================
* 3. MERGE ALL VARIABLES INTO ONE PANEL
* ==============================================================================
* Each of the three .dta files saved above has one row per region-year pair
* and holds exactly one variable (plus region_id and year). This section
* combines them into a single panel with one row per region-year pair and
* all three variables as columns.
*
* personnel.dta is loaded first as the starting point. Each subsequent
* merge is a 1:1 match on region_id and year: for every region-year row in
* the file currently in memory, Stata looks up the matching region-year row
* in the "using" file and attaches its variable. nogenerate suppresses the
* _merge indicator variable that Stata creates by default, since with only
* three files (all built from the same 20 regions x 46 years grid) a full
* match is expected every time.

local allvars personnel out_migr fert_rate

use "${reshaped_dir}/personnel.dta", clear

foreach v of local allvars {
    if "`v'" != "personnel" {
        merge 1:1 region_id year using "${reshaped_dir}/`v'.dta", nogenerate
    }
}

* ------------------------------------------------------------------------------
* 3.1: sort and declare the panel structure
* ------------------------------------------------------------------------------
* xtset tells Stata that region_id identifies the cross-sectional unit
* (each region) and year identifies the time dimension. This is required
* before using panel-data commands such as ipolate (section 3.2) or synth.

sort  region_id year
xtset region_id year

* ------------------------------------------------------------------------------
* 3.2: fix the known gap in personnel (2014-2015)
* ------------------------------------------------------------------------------
* ISTAT indicator 9110 (personnel) is missing for the years 2014 and 2015
* for every one of the 20 regions. This is a gap in ISTAT's own reporting, not
* a region-specific data problem, so it is reasonable to fill it in rather
* than lose two years of the panel.
*
* ipolate fills each missing year with a straight line drawn between the
* nearest observed value before it (2013) and the nearest observed value
* after it (2016). In other words:
*   personnel(2014) = personnel(2013) + 1/3 x [personnel(2016) - personnel(2013)]
*   personnel(2015) = personnel(2013) + 2/3 x [personnel(2016) - personnel(2013)]
* "bysort region_id" makes ipolate do this separately for each region.
*
* The replace is limited to years 2014-2015 AND to cells that are still
* missing, so that if any region already has a real value in those years,
* it is not overwritten by the interpolated one.

bysort region_id (year): ipolate personnel year, gen(_temp_pers)
replace personnel = _temp_pers if inlist(year, 2014, 2015) & missing(personnel)
drop _temp_pers

* These "assert" lines are safety checks: if any statement is false, Stata
* stops immediately with an error instead of silently continuing with bad
* data. This confirms the interpolation above worked as intended and that
* out_migr and fert_rate genuinely have no gaps within their available
* windows.

assert !missing(personnel) if inrange(year, 1994, 2022)
assert !missing(out_migr)  if inrange(year, 1999, 2023)
assert !missing(fert_rate) if inrange(year, 1980, 2024)



* ==============================================================================
* 4. UNIT HARMONISATION
* ==============================================================================

* ------------------------------------------------------------------------------
* 4.1: personnel -- rescale from per 10,000 to per 1,000 inhabitants
* ------------------------------------------------------------------------------
* ISTAT 9110 reports personnel per 10,000 inhabitants. Dividing by 10 puts
* it on a per-1,000 basis.

replace personnel = personnel / 10

* ------------------------------------------------------------------------------
* 4.2: fert_rate -- rescale from per 1,000 women to per woman
* ------------------------------------------------------------------------------
* ISTAT 0180 reports the total fertility rate as children per 1,000 women.
* Dividing by 1,000 converts it to the standard "children per woman" scale,
* which is how this indicator is conventionally reported in the literature.

replace fert_rate = fert_rate / 1000



* ==============================================================================
* 5. TREATMENT AND DONOR-POOL INDICATORS
* ==============================================================================

* ------------------------------------------------------------------------------
* 5.1: active PdR dummy (region x year)
* ------------------------------------------------------------------------------
* pdr = 1 marks a region-year in which that region's recovery plan was active. 
* Built by hand from each region's known entry and (where applicable) exit year.

gen pdr = 0

* Lazio, Abruzzo, Molise, Campania, Sicilia: PdR from 2007 onward
replace pdr = 1 if inlist(region_id, 12, 13, 14, 15, 19) & year >= 2007

* Liguria and Sardegna: PdR 2007-2009 only (exit after 2009)
replace pdr = 1 if inlist(region_id, 7, 20) & inrange(year, 2007, 2009)

* Calabria: PdR from 2009 onward
replace pdr = 1 if region_id == 18 & year >= 2009

* Piemonte: PdR 2010-2016 (exits 2017)
replace pdr = 1 if region_id ==  1 & inrange(year, 2010, 2016)

* Puglia: PdR from 2010 onward
replace pdr = 1 if region_id == 16 & year >= 2010

label var pdr "=1 if region-year is under an active PdR"

* ------------------------------------------------------------------------------
* 5.2: donor-pool dummy
* ------------------------------------------------------------------------------
* donor = 1 marks the ten regions that were never placed under a PdR at any
* point in the sample. These are the only regions eligible to be used as
* SCM donors.
*
* Never-treated regions: Valle d'Aosta (2), Lombardia (3), Trentino-A.A. (4),
* Veneto (5), Friuli-V.G. (6), Emilia-Romagna (8), Toscana (9), Umbria (10),
* Marche (11), Basilicata (17).

gen donor = inlist(region_id, 2, 3, 4, 5, 6, 8, 9, 10, 11, 17)
label var donor "=1 if region is in the SCM donor pool (never treated)"



* ==============================================================================
* 6. VARIABLE LABELS AND FINAL ORDERING
* ==============================================================================

label var region_id "ISTAT region code (1-20)"
label var year      "Calendar year"

label var personnel "SSN personnel rate, per 1,000 inhabitants (ISTAT 9110)"
label var out_migr  "Patient out-migration rate, % of discharged residents (ISTAT 7343)"
label var fert_rate "Total fertility rate, TFT (children per woman) (ISTAT 0180) -- placebo"

* Reorder columns: identifiers and treatment indicators first, then the two
* outcomes, then the placebo variable.
order region_id year pdr donor personnel out_migr fert_rate



* ==============================================================================
* 7. VERIFY PANEL STRUCTURE AND SAVE
* ==============================================================================
* xtdescribe summarises the panel's shape (regions, years, balance) so any
* unexpected gap or duplicate would show up here before saving. The two
* tabulate commands are spot-checks confirming that the pdr indicator is
* coded as intended: exactly the expected regions and years show pdr = 1.

xtdescribe

tabulate year      pdr, missing
tabulate region_id pdr, missing

save "${out_dir}/panel_HFA.dta", replace



* ==============================================================================
* END OF .DO FILE
* ==============================================================================

log close
