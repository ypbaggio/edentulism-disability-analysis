version 19.0
clear all
set maxvar 32767
set more off
set linesize 200
capture log close _all

/*
Participant-level chained multiple imputation for the six-cohort primary
two-stage analysis. The completed participant records are merged back into the
first-event interval panel only inside this local Stata session. Primary outcome
models are unweighted; observation weighting is retained as a separate
missingness-indicator sensitivity analysis. Only aggregate cohort estimates are
exported.

The public workflow begins with a harmonized, analysis-ready first-event panel.
Raw cohort files and source-data construction are outside this repository.

Usage:
  do 01_primary_mi_models.do "harmonized_first_event_panel.dta" "output_directory"
*/

args input_file output_dir
if `"`input_file'"' == "" | `"`output_dir'"' == "" {
    display as error "Supply the combined first-event panel and output directory."
    exit 198
}

capture mkdir `"`output_dir'"'
log using `"`output_dir'/JDR_R2_multiple_imputation.log"', text replace
use `"`input_file'"', clear

foreach variable in participant_id cohort interval_number analysis_time first_disability ///
    death_before_disability edentulous age_start female education_raw married_raw ///
    current_smoker_raw comorbidity_count_raw {
    confirm variable `variable'
}
assert !missing(participant_id, cohort, interval_number, analysis_time, edentulous, age_start, female)
assert analysis_time > 0
isid participant_id interval_number

tempfile interval_panel mi_baseline mi_fits
save `interval_panel'

bysort participant_id: egen double total_followup = total(analysis_time)
bysort participant_id: egen byte disability_ever = max(first_disability)
bysort participant_id: egen byte death_ever = max(death_before_disability)
bysort participant_id (interval_number): keep if _n == 1
encode cohort, gen(cohort_id)

preserve
gen byte one_person = 1
foreach variable in education_raw married_raw current_smoker_raw comorbidity_count_raw {
    gen byte `variable'_missing = missing(`variable')
}
collapse (sum) n_people=one_person education_missing=education_raw_missing ///
    married_missing=married_raw_missing smoker_missing=current_smoker_raw_missing ///
    comorbidity_missing=comorbidity_count_raw_missing, by(cohort)
foreach variable in education married smoker comorbidity {
    gen double `variable'_missing_percent = 100 * `variable'_missing / n_people
}
export delimited using `"`output_dir'/JDR_multiple_imputation_missingness.csv"', replace
restore

keep participant_id cohort cohort_id edentulous age_start female education_raw ///
    married_raw current_smoker_raw comorbidity_count_raw total_followup ///
    disability_ever death_ever

mi set mlong
mi register imputed education_raw married_raw current_smoker_raw comorbidity_count_raw
mi register regular participant_id cohort cohort_id edentulous age_start female ///
    total_followup disability_ever death_ever

mi impute chained ///
    (pmm, knn(5)) education_raw comorbidity_count_raw ///
    (logit) married_raw current_smoker_raw ///
    = edentulous c.age_start i.female c.total_followup ///
      i.disability_ever i.death_ever, ///
    add(20) burnin(10) rseed(20260927) by(cohort_id) augment dots

mi describe
save `mi_baseline'

forvalues m = 1/20 {
    preserve
    mi extract `m', clear
    keep participant_id education_raw married_raw current_smoker_raw comorbidity_count_raw
    isid participant_id
    save `"`output_dir'/LOCAL_ONLY_MI_baseline_`m'.dta"', replace
    restore
}

use `interval_panel', clear
levelsof cohort, local(cohorts)
tempname fits
postfile `fits' str8 cohort str28 outcome int imputation double log_irr within_variance ///
    n_people n_events using `mi_fits', replace

forvalues m = 1/20 {
    use `interval_panel', clear
    drop education_raw married_raw current_smoker_raw comorbidity_count_raw
    merge m:1 participant_id using `"`output_dir'/LOCAL_ONLY_MI_baseline_`m'.dta"', ///
        assert(match) keep(match) nogen
    egen long person_numeric = group(participant_id)

    foreach current_cohort of local cohorts {
        preserve
        keep if cohort == "`current_cohort'"
        quietly summarize education_raw
        if r(sd) <= 0 | missing(r(sd)) {
            display as error "Education has no variation in `current_cohort', imputation `m'."
            exit 459
        }
        gen double education_z_mi = (education_raw - r(mean)) / r(sd)
        quietly summarize comorbidity_count_raw
        if r(sd) <= 0 | missing(r(sd)) {
            display as error "Comorbidity count has no variation in `current_cohort', imputation `m'."
            exit 459
        }
        gen double comorbidity_z_mi = (comorbidity_count_raw - r(mean)) / r(sd)
        mkspline age_rcs_mi = age_start, cubic nknots(4)
        foreach outcome in first_disability death_before_disability first_exit {
            capture quietly poisson `outcome' edentulous age_rcs_mi* female ///
                education_z_mi married_raw current_smoker_raw comorbidity_z_mi ///
                i.interval_number, exposure(analysis_time) ///
                vce(cluster person_numeric) difficult iterate(200)
            if _rc {
                local poisson_rc = _rc
                display as text "Poisson ML did not converge for `current_cohort', imputation `m', `outcome' (r(`poisson_rc')); retrying with IRLS."
                capture quietly glm `outcome' edentulous age_rcs_mi* female ///
                    education_z_mi married_raw current_smoker_raw comorbidity_z_mi ///
                    i.interval_number, family(poisson) link(log) ///
                    exposure(analysis_time) vce(cluster person_numeric) ///
                    iterate(200)
                if _rc {
                    local glm_rc = _rc
                    display as error "Both Poisson estimators failed for `current_cohort', imputation `m', `outcome' (poisson r(`poisson_rc'); glm r(`glm_rc'))."
                    exit `glm_rc'
                }
            }
            scalar beta = _b[edentulous]
            scalar beta_variance = _se[edentulous]^2
            quietly egen byte __tag = tag(person_numeric) if e(sample)
            quietly count if __tag == 1
            scalar n_person = r(N)
            drop __tag
            quietly summarize `outcome' if e(sample), meanonly
            scalar n_event = r(sum)
            post `fits' ("`current_cohort'") ("`outcome'") (`m') (beta) ///
                (beta_variance) (n_person) (n_event)
        }
        restore
    }
}
postclose `fits'

use `mi_fits', clear
isid cohort outcome imputation
bysort cohort outcome: egen int m_used = count(log_irr)
bysort cohort outcome: egen double pooled_log_irr = mean(log_irr)
bysort cohort outcome: egen double within_mean = mean(within_variance)
bysort cohort outcome: egen double between_sd = sd(log_irr)
gen double between_variance = between_sd^2
gen double total_variance = within_mean + (1 + 1 / m_used) * between_variance
gen double se = sqrt(total_variance)
gen double rubin_df = (m_used - 1) * ///
    (1 + within_mean / ((1 + 1 / m_used) * between_variance))^2
replace rubin_df = 1000000 if between_variance < 1e-12 | missing(rubin_df)
gen double critical = invttail(rubin_df, .025)
gen double irr = exp(pooled_log_irr)
gen double lcl = exp(pooled_log_irr - critical * se)
gen double ucl = exp(pooled_log_irr + critical * se)
gen double p_value = 2 * ttail(rubin_df, abs(pooled_log_irr / se))
bysort cohort outcome (imputation): keep if _n == 1
rename pooled_log_irr log_irr_mi
keep cohort outcome m_used log_irr_mi se irr lcl ucl p_value rubin_df ///
    within_mean between_variance total_variance n_people n_events
order cohort outcome m_used log_irr_mi se irr lcl ucl p_value rubin_df ///
    n_people n_events within_mean between_variance total_variance
format log_irr_mi se within_mean between_variance total_variance %10.5f
format irr lcl ucl %9.3f
format p_value %9.4f
export delimited using `"`output_dir'/JDR_multiple_imputation_cohort_results.csv"', replace
list, noobs clean

forvalues m = 1/20 {
    erase `"`output_dir'/LOCAL_ONLY_MI_baseline_`m'.dta"'
}

display as result "Multiple-imputation cohort models completed."
display as text "Only aggregate results were exported; participant-level imputations remain temporary."
log close
