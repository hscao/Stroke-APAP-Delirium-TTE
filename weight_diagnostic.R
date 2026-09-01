library(survival)
library(survivalCCW)
library(dplyr)
library(openxlsx)

data <- read.xlsx("apap_delirium.xlsx")

data <- data %>%
    mutate(
        exposure = as.integer(exposure),
        event = as.integer(event)
    )

covars <- c(
    "age", "sex", "aud", "hypertension", "diabetes", "hyperlipidemia",
    "previous_stroke", "chf", "copd", "ckd", "dementia", "aphasia",
    "infection", "opioids", "benzodiazepines", "dexmedetomidine",
    "glucocorticoid", "mv", "gcs", "stroke_type1", "stroke_type2"
)

clones <- create_clones(
    df = data,
    id = "subject_id",
    event = "event",
    time_to_event = "timetoevent",
    exposure = "exposure",
    time_to_exposure = "timetoexposure",
    ced_window = 48
)

clones_long <- cast_clones_to_long(clones)

untruncated_weights <- generate_ccw(clones_long, predvars = covars)

min(untruncated_weights$weight_cox)
# 1
median(untruncated_weights$weight_cox)
# 1
max(untruncated_weights$weight_cox)
# 55.53241


# Winsorization thresholds
thresholds <- c(Untruncated = 1, P99 = 0.99, P97.5 = 0.975, P95 = 0.95, P90 = 0.90, P85 = 0.85)

kish_ess <- function(w) sum(w)^2 / sum(w^2)

weight_diagnostics <- do.call(rbind, lapply(names(thresholds), function(label) {
    weights <- if (label == "Untruncated") {
        untruncated_weights
    } else {
        winsorize_ccw_weights(untruncated_weights, quantiles = c(0, thresholds[[label]]))
    }

    terminal_weights <- weights %>%
        arrange(subject_id, clone, t_stop) %>%
        group_by(subject_id, clone) %>%
        slice_tail(n = 1) %>%
        ungroup()

    data.frame(
        `Weight handling` = label,
        `Weights modified` = 100 * mean(abs(weights$weight_cox - untruncated_weights$weight_cox) > 1e-12),
        `Maximum weight` = max(weights$weight_cox),
        ESS = kish_ess(terminal_weights$weight_cox),
        check.names = FALSE
    )
}))

print(weight_diagnostics, row.names = FALSE)
#  Weight handling Weights modified Maximum weight      ESS
#      Untruncated         0.000000      55.532410 2305.160
#              P99         1.011170      21.934059 2726.315
#            P97.5         2.499838      15.645361 3045.796
#              P95         4.981423      10.120365 4039.654
#              P90         9.944290       3.308498 5811.977
#              P85        14.609545       2.315362 6283.419
