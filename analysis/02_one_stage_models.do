version 19.0
clear all
set maxvar 32767
set more off
set linesize 200

/*
Secondary one-stage interval-level IPD analysis.

Usage:
  do 02_one_stage_models.do "harmonized_first_event_panel.dta" "output_directory"

The input panel must contain one row per observed risk interval and the
variables listed in data/input_schema.csv. It is constructed locally under
the applicable cohort access and data-use terms.
*/

args input_file output_dir
if `"`input_file'"' == "" {
    display as error "Supply the harmonized first-event panel path as argument 1."
    exit 198
}
if `"`output_dir'"' == "" local output_dir "."

capture mkdir `"`output_dir'"'
log using `"`output_dir'/JDR_R2_one_stage_model.log"', text replace
use `"`input_file'"', clear

assert !missing(participant_id, cohort, analysis_time, edentulous, age_start)
assert analysis_time > 0
assert inlist(edentulous, 0, 1)
assert inlist(first_disability, 0, 1)
assert inlist(death_before_disability, 0, 1)
assert inlist(first_exit, 0, 1)

replace share_country = "Not SHARE" if cohort != "SHARE"
replace analysis_unit = cohort if cohort != "SHARE"
encode cohort, gen(cohort_id)
egen long person_numeric = group(participant_id)
quietly tabulate share_country if cohort == "SHARE", generate(share_ctry_)
foreach variable of varlist share_ctry_* {
    replace `variable' = 0 if missing(`variable')
}
drop share_ctry_1

capture drop age_rcs*
mkspline age_rcs = age_start, cubic nknots(4)

local cohort_age_terms
foreach agevar of varlist age_rcs* {
    local cohort_age_terms `cohort_age_terms' c.`agevar'#i.cohort_id
}
local covars edentulous `cohort_age_terms' female education_z married_imp smoker_imp comorbidity_z education_missing married_missing smoker_missing comorbidity_missing i.interval_number share_ctry_*
scalar tcrit5 = invttail(5, .025)

tempname post_handle
postfile `post_handle' str28 outcome str48 model str32 inference double log_irr se irr lcl ucl p_value wald_lcl wald_ucl n_people n_events using `"`output_dir'/JDR_one_stage_results.dta"', replace

foreach outcome in first_disability death_before_disability first_exit {
    capture mepoisson `outcome' `covars', exposure(analysis_time) ///
        || cohort_id:, intpoints(7) difficult iterate(200)
    local fit_rc = _rc
    local retry_fit = (`fit_rc' != 0)
    if `fit_rc' == 0 {
        if e(converged) != 1 local retry_fit = 1
    }
    if `retry_fit' {
        display as text "The initial one-stage `outcome' fit did not converge; retrying with 12 integration points."
        capture mepoisson `outcome' `covars', exposure(analysis_time) ///
            || cohort_id:, intpoints(12) difficult iterate(500)
        local fit_rc = _rc
    }
    local failed_fit = (`fit_rc' != 0)
    if `fit_rc' == 0 {
        if e(converged) != 1 local failed_fit = 1
    }
    if `failed_fit' {
        display as error "The one-stage `outcome' model did not converge after the planned retry."
        exit 430
    }
    scalar beta = _b[edentulous]
    scalar se = _se[edentulous]
    tempvar sample_person
    quietly egen byte `sample_person' = tag(participant_id) if e(sample)
    quietly count if `sample_person' == 1
    scalar n_person = r(N)
    drop `sample_person'
    quietly summarize `outcome' if e(sample), meanonly
    scalar n_event = r(sum)
    post `post_handle' ("`outcome'") ("cohort_RE_country_FE_cohort_age") ("t distribution, 5 df") (beta) (se) ///
        (exp(beta)) (exp(beta-tcrit5*se)) (exp(beta+tcrit5*se)) ///
        (2*ttail(5, abs(beta/se))) (exp(beta-invnormal(.975)*se)) (exp(beta+invnormal(.975)*se)) ///
        (n_person) (n_event)
    estimates store one_`outcome'
}
postclose `post_handle'

preserve
use `"`output_dir'/JDR_one_stage_results.dta"', clear
format log_irr se %9.4f
format irr lcl ucl %9.3f
format p_value %9.4f
export delimited using `"`output_dir'/JDR_one_stage_results.csv"', replace
list, noobs clean
restore

display as text "One-stage models completed. Review convergence and variance components in the log."
log close
