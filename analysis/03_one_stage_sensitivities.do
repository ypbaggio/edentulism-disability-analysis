version 19.0
clear all
set maxvar 32767
set more off
set linesize 200

/*
One-stage sensitivity analyses.

Usage:
  do 03_one_stage_sensitivities.do "harmonized_first_event_panel.dta" "output_directory"

This script complements, but does not replace, the primary mixed-effects
one-stage model. It estimates cohort-fixed and SHARE-country-fixed models with
participant-clustered covariance, first without observation weights and then
with the cohort-specific inverse-probability-of-observation weights.
*/

args input_file output_dir
if `"`input_file'"' == "" {
    display as error "Supply the harmonized first-event panel path as argument 1."
    exit 198
}
if `"`output_dir'"' == "" local output_dir "."

capture mkdir `"`output_dir'"'
capture log close _all
log using `"`output_dir'/JDR_R2_one_stage_sensitivity.log"', text replace
use `"`input_file'"', clear

assert !missing(participant_id, cohort, analysis_time, edentulous, age_start, ipow)
assert analysis_time > 0
assert ipow > 0
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
local covars edentulous `cohort_age_terms' female education_z married_imp smoker_imp comorbidity_z education_missing married_missing smoker_missing comorbidity_missing i.interval_number i.cohort_id share_ctry_*
local mixed_covars edentulous `cohort_age_terms' female education_z married_imp smoker_imp comorbidity_z education_missing married_missing smoker_missing comorbidity_missing i.interval_number share_ctry_*
scalar tcrit5 = invttail(5, .025)

tempname post_handle
postfile `post_handle' str28 outcome str48 model double log_irr se irr lcl ucl p_value n_people n_events n_intervals using `"`output_dir'/JDR_one_stage_sensitivity_results.dta"', replace

foreach outcome in first_disability death_before_disability first_exit {
    capture quietly poisson `outcome' `covars', exposure(analysis_time) ///
        vce(cluster person_numeric) difficult iterate(200)
    if _rc {
        local poisson_rc = _rc
        display as text "Unweighted Poisson ML did not converge for `outcome' (r(`poisson_rc')); retrying with IRLS."
        capture quietly glm `outcome' `covars', family(poisson) link(log) ///
            exposure(analysis_time) vce(cluster person_numeric) iterate(200)
        if _rc {
            local glm_rc = _rc
            display as error "Both unweighted estimators failed for `outcome' (poisson r(`poisson_rc'); glm r(`glm_rc'))."
            exit `glm_rc'
        }
    }
    scalar beta = _b[edentulous]
    scalar stderr = _se[edentulous]
    tempvar sample_person
    quietly egen byte `sample_person' = tag(participant_id) if e(sample)
    quietly count if `sample_person' == 1
    scalar n_person = r(N)
    drop `sample_person'
    quietly summarize `outcome' if e(sample), meanonly
    scalar n_event = r(sum)
    scalar n_interval = e(N)
    post `post_handle' ("`outcome'") ("cohort_country_FE_clustered") (beta) (stderr) ///
        (exp(beta)) (exp(beta-invnormal(.975)*stderr)) (exp(beta+invnormal(.975)*stderr)) ///
        (2*normal(-abs(beta/stderr))) (n_person) (n_event) (n_interval)

    capture quietly poisson `outcome' `covars' [pw=ipow], exposure(analysis_time) ///
        vce(cluster person_numeric) difficult iterate(200)
    if _rc {
        local poisson_rc = _rc
        display as text "IPOW Poisson ML did not converge for `outcome' (r(`poisson_rc')); retrying with IRLS."
        capture quietly glm `outcome' `covars' [pw=ipow], family(poisson) link(log) ///
            exposure(analysis_time) vce(cluster person_numeric) iterate(200)
        if _rc {
            local glm_rc = _rc
            display as error "Both IPOW estimators failed for `outcome' (poisson r(`poisson_rc'); glm r(`glm_rc'))."
            exit `glm_rc'
        }
    }
    scalar beta = _b[edentulous]
    scalar stderr = _se[edentulous]
    tempvar sample_person
    quietly egen byte `sample_person' = tag(participant_id) if e(sample)
    quietly count if `sample_person' == 1
    scalar n_person = r(N)
    drop `sample_person'
    quietly summarize `outcome' if e(sample), meanonly
    scalar n_event = r(sum)
    scalar n_interval = e(N)
    post `post_handle' ("`outcome'") ("cohort_country_FE_IPOW_clustered") (beta) (stderr) ///
        (exp(beta)) (exp(beta-invnormal(.975)*stderr)) (exp(beta+invnormal(.975)*stderr)) ///
        (2*normal(-abs(beta/stderr))) (n_person) (n_event) (n_interval)
}

capture quietly poisson first_disability `covars' [pw=ipow] if ///
    education_missing == 0 & married_missing == 0 & smoker_missing == 0 & ///
    comorbidity_missing == 0, exposure(analysis_time) ///
    vce(cluster person_numeric) difficult iterate(200)
if _rc {
    local poisson_rc = _rc
    display as text "Complete-case Poisson ML did not converge (r(`poisson_rc')); retrying with IRLS."
    capture quietly glm first_disability `covars' [pw=ipow] if ///
        education_missing == 0 & married_missing == 0 & smoker_missing == 0 & ///
        comorbidity_missing == 0, family(poisson) link(log) ///
        exposure(analysis_time) vce(cluster person_numeric) iterate(200)
    if _rc {
        local glm_rc = _rc
        display as error "Both complete-case estimators failed (poisson r(`poisson_rc'); glm r(`glm_rc'))."
        exit `glm_rc'
    }
}
scalar beta = _b[edentulous]
scalar stderr = _se[edentulous]
tempvar sample_person
quietly egen byte `sample_person' = tag(participant_id) if e(sample)
quietly count if `sample_person' == 1
scalar n_person = r(N)
drop `sample_person'
quietly summarize first_disability if e(sample), meanonly
scalar n_event = r(sum)
scalar n_interval = e(N)
post `post_handle' ("first_disability") ("complete_case_IPOW_clustered") (beta) (stderr) ///
    (exp(beta)) (exp(beta-invnormal(.975)*stderr)) (exp(beta+invnormal(.975)*stderr)) ///
    (2*normal(-abs(beta/stderr))) (n_person) (n_event) (n_interval)

capture quietly mepoisson first_disability `mixed_covars', exposure(analysis_time) ///
    || cohort_id: edentulous, covariance(unstructured) intpoints(7) difficult
if _rc == 0 & e(converged) == 1 {
    scalar beta = _b[edentulous]
    scalar stderr = _se[edentulous]
    tempvar sample_person
    quietly egen byte `sample_person' = tag(participant_id) if e(sample)
    quietly count if `sample_person' == 1
    scalar n_person = r(N)
    drop `sample_person'
    quietly summarize first_disability if e(sample), meanonly
    scalar n_event = r(sum)
    scalar n_interval = e(N)
    post `post_handle' ("first_disability") ("cohort_RE_random_slope_t5") (beta) (stderr) ///
        (exp(beta)) (exp(beta-tcrit5*stderr)) (exp(beta+tcrit5*stderr)) ///
        (2*ttail(5, abs(beta/stderr))) (n_person) (n_event) (n_interval)
}
else {
    display as text "Random-slope sensitivity did not converge and was not exported."
}
postclose `post_handle'

preserve
use `"`output_dir'/JDR_one_stage_sensitivity_results.dta"', clear
format log_irr se %9.4f
format irr lcl ucl %9.3f
format p_value %9.4f
format n_people n_events n_intervals %12.0f
export delimited using `"`output_dir'/JDR_one_stage_sensitivity_results.csv"', replace
list, noobs clean
restore

display as text "One-stage fixed-effects sensitivity models completed."
log close
