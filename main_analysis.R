library(survival)
library(survivalCCW)
library(tidyverse)
library(openxlsx)
library(survminer)
library(boot)
library(forestplot)
library(EValue)


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



# bootstrap
boot_cox <- function(data, indices) {
    ccw_df <- data[indices, ] |>
        create_clones(
            id = "subject_id",
            event = "event",
            time_to_event = "timetoevent",
            exposure = "exposure",
            time_to_exposure = "timetoexposure",
            ced_window = 48
        ) |>
        cast_clones_to_long() |>
        generate_ccw(predvars = covars) |>
        winsorize_ccw_weights(quantiles = c(0, 0.90))

    cox_ccw <- coxph(Surv(t_start, t_stop, outcome) ~ clone,
        data = ccw_df, weights = weight_cox
    )

    # HR
    hr <- cox_ccw |>
        coef() |>
        exp()
    out <- c("hr" = hr)

    surv_1 <- survfit(Surv(t_start, t_stop, outcome) ~ 1L,
        data = ccw_df[ccw_df$clone == 1, ],
        weights = weight_cox
    )
    surv_0 <- survfit(Surv(t_start, t_stop, outcome) ~ 1L,
        data = ccw_df[ccw_df$clone == 0, ],
        weights = weight_cox
    )

    # RMST & RMTL
    rmst_1 <- surv_1 %>%
        summary(rmean = 120) %>%
        pluck("table") %>%
        pluck("rmean")
    rmst_0 <- surv_0 %>%
        summary(rmean = 120) %>%
        pluck("table") %>%
        pluck("rmean")

    rmtl_1 <- 120 - rmst_1
    rmtl_0 <- 120 - rmst_0

    rmtl_diff <- rmtl_1 - rmtl_0
    out <- c(out, "rmtl_diff" = rmtl_diff)
    out

    # surv_diff & inc_diff
    surv_1_120hr <- summary(surv_1, times = 120, extend = TRUE)$surv[[1]]
    surv_0_120hr <- summary(surv_0, times = 120, extend = TRUE)$surv[[1]]

    inc_1_120hr <- 1 - surv_1_120hr
    inc_0_120hr <- 1 - surv_0_120hr

    inc_diff_120hr <- inc_1_120hr - inc_0_120hr
    out <- c(out, "inc_diff_120hr" = inc_diff_120hr)
}

# 500 bootstrap resampling
R_main <- 500
pb_main <- txtProgressBar(min = 0, max = R_main, style = 3)
iter_main <- 0L
stat_with_progress <- function(d, idx) {
    n <- NROW(d)
    is_t0 <- length(idx) == n && identical(idx, seq_len(n))
    if (!is_t0) {
        iter_main <<- iter_main + 1L
        setTxtProgressBar(pb_main, iter_main, label = sprintf("[Main] %d/%d", iter_main, R_main))
    }
    boot_cox(d, idx)
}
start_main <- Sys.time()
boot_out <- boot(data = data, statistic = stat_with_progress, R = R_main)
close(pb_main)
cat("\n")

elapsed_main <- as.numeric(difftime(Sys.time(), start_main, units = "mins"))
cat(sprintf("[Main] Complete | Total time %.1f minutes\n", elapsed_main))
flush.console()

hr_ci <- boot.ci(boot_out, type = "norm", index = 1)
rmtl_ci <- boot.ci(boot_out, type = "norm", index = 2)
inc_ci <- boot.ci(boot_out, type = "norm", index = 3)

hr_est <- boot_out$t0[1]
rmtl_est <- boot_out$t0[2]
inc_est <- boot_out$t0[3]

cat(sprintf(
    "HR = %.4f (95%% CI: %.4f, %.4f)\n",
    hr_est, hr_ci$normal[2], hr_ci$normal[3]
))
cat(sprintf(
    "RMTL difference = %.4f (95%% CI: %.4f, %.4f)\n",
    rmtl_est, rmtl_ci$normal[2], rmtl_ci$normal[3]
))
cat(sprintf(
    "Risk difference at 120h = %.4f (95%% CI: %.4f, %.4f)\n",
    inc_est, inc_ci$normal[2], inc_ci$normal[3]
))



# Weighted cumulative incidence curve
weights_winsorized <- data |>
    create_clones(
        id = "subject_id",
        event = "event",
        time_to_event = "timetoevent",
        exposure = "exposure",
        time_to_exposure = "timetoexposure",
        ced_window = 48
    ) |>
    cast_clones_to_long() |>
    generate_ccw(predvars = covars) |>
    winsorize_ccw_weights(quantiles = c(0, 0.90))

weighted_fit <- survfit(Surv(t_start, t_stop, outcome) ~ clone, data = weights_winsorized, weights = weight_cox)

km_plot <- ggsurvplot(weighted_fit,
    data = weights_winsorized,
    fun = "event",
    linetype = 1,
    size = 1.2,
    conf.int = FALSE,
    conf.int.alpha = 0.3,
    xlab = "Time of follow-up, h",
    ylab = "Cumulative 5-day delirium incidence, %",
    legend.title = "",
    legend.labs = c("No acetaminophen group", "Acetaminophen group"),
    legend = c(0.30, 0.80),
    break.x.by = 24,
    xlim = c(0, 122),
    break.y.by = 0.1,
    ylim = c(0.0, 1.0),
    axes.offset = FALSE,
    surv.median.line = "none",
    ggtheme = theme_survminer() +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            panel.grid.major.y = element_line(color = "#e2e2e3", linetype = "solid"),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16, hjust = 0.5)
        ),
    table.theme = theme_bw(),
    palette = c("#ff9027", "#2b5461"),
    pval = FALSE,
    pval.method = FALSE,
    censor = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.22,
    fontsize = 5,
    tables.y.text = FALSE,
    tables.y.text.col = TRUE
)

km_plot$plot <- km_plot$plot +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.1),
        labels = scales::number_format(scale = 100, accuracy = 1),
        expand = c(0, 0)
    )

km_plot



# Subgroup analysis
boot_sub <- function(data, indices) {
    sample_df <- data[indices, ]

    sub_df <- tryCatch(
        {
            sample_df |>
                create_clones(
                    id = "subject_id",
                    event = "event",
                    time_to_event = "timetoevent",
                    exposure = "exposure",
                    time_to_exposure = "timetoexposure",
                    ced_window = 48
                ) |>
                cast_clones_to_long() |>
                generate_ccw(predvars = covars) |>
                winsorize_ccw_weights(quantiles = c(0, 0.90))
        },
        error = function(e) NULL
    )

    if (is.null(sub_df) || nrow(sub_df) == 0) {
        return(NA_real_)
    }

    events_by_arm <- tapply(sub_df$outcome, sub_df$clone, function(x) sum(x == 1, na.rm = TRUE))
    if (length(events_by_arm) < 2 || any(is.na(events_by_arm)) || any(events_by_arm == 0)) {
        return(NA_real_)
    }

    fit <- suppressWarnings(tryCatch(
        coxph(Surv(t_start, t_stop, outcome) ~ clone,
            data = sub_df, weights = weight_cox
        ),
        error = function(e) NULL
    ))

    if (is.null(fit)) {
        return(NA_real_)
    }

    beta <- coef(fit)[["clone"]]

    if (!is.finite(beta)) {
        return(NA_real_)
    }

    beta
}

subgroup_list <- list(
    "<65" = subset(data, age < 65),
    "≥65" = subset(data, age >= 65),
    "Male" = subset(data, sex == 1),
    "Female" = subset(data, sex == 0),
    "Ischemic stroke" = subset(data, stroke_type == 1),
    "ICH" = subset(data, stroke_type == 2),
    "SAH" = subset(data, stroke_type == 3),
    "Yes CKD" = subset(data, ckd == 1),
    "No CKD" = subset(data, ckd == 0)
)

all_groups_list <- subgroup_list
bootstrap_results <- list()
R_sub <- 500
global_start <- Sys.time()

make_statistic <- function(subgroup_name, R_sub, pb) {
    iter <- 0L
    start_time <- Sys.time()
    function(d, idx) {
        is_t0 <- length(idx) == NROW(d) && identical(idx, seq_len(NROW(d)))
        if (!is_t0) {
            iter <<- iter + 1L
            elapsed_sub <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
            avg_sub <- elapsed_sub / iter
            eta_sub <- max(0, avg_sub * (R_sub - iter))
            setTxtProgressBar(pb, iter, label = sprintf("[%s] %d/%d | ETA %.1fm", subgroup_name, iter, R_sub, eta_sub))
        }
        boot_sub(d, idx)
    }
}

for (subgroup_name in names(all_groups_list)) {
    cat(paste("Running bootstrap for subgroup:", subgroup_name, "...\n"))
    subgroup_data <- all_groups_list[[subgroup_name]]

    subgroup_start <- Sys.time()
    pb <- txtProgressBar(min = 0, max = R_sub, style = 3)
    stat_fn <- make_statistic(subgroup_name, R_sub, pb)

    boot_out_sub <- boot(
        data = subgroup_data,
        statistic = stat_fn,
        R = R_sub
    )

    close(pb)
    cat("\n")
    cat(sprintf(
        "[%s] Complete | Total time %.1fm\n",
        subgroup_name,
        as.numeric(difftime(Sys.time(), subgroup_start, units = "mins"))
    ))
    flush.console()

    point_estimate <- boot_out_sub$t0[1]
    ci <- boot.ci(boot_out_sub, type = "norm", index = 1)

    hr_point <- exp(point_estimate)
    hr_lower <- exp(ci$normal[2])
    hr_upper <- exp(ci$normal[3])

    bootstrap_results[[subgroup_name]] <- list(
        hr = hr_point,
        lower = hr_lower,
        upper = hr_upper
    )
}

cat(sprintf(
    "All subgroups complete | Total time %.1fm\n",
    as.numeric(difftime(Sys.time(), global_start, units = "mins"))
))
flush.console()

n_vec <- sapply(all_groups_list, nrow)

percent_vec <- c(
    round(n_vec["<65"] / sum(!is.na(data$age)) * 100, 1),
    round(n_vec["≥65"] / sum(!is.na(data$age)) * 100, 1),
    round(n_vec["Male"] / sum(!is.na(data$sex)) * 100, 1),
    round(n_vec["Female"] / sum(!is.na(data$sex)) * 100, 1),
    round(n_vec["Ischemic stroke"] / sum(!is.na(data$stroke_type)) * 100, 1),
    round(n_vec["ICH"] / sum(!is.na(data$stroke_type)) * 100, 1),
    round(n_vec["SAH"] / sum(!is.na(data$stroke_type)) * 100, 1),
    round(n_vec["Yes CKD"] / sum(!is.na(data$ckd)) * 100, 1),
    round(n_vec["No CKD"] / sum(!is.na(data$ckd)) * 100, 1)
)

results_df <- do.call(rbind, lapply(names(bootstrap_results), function(name) {
    data.frame(
        subgroup = name,
        hr = bootstrap_results[[name]]$hr,
        lower = bootstrap_results[[name]]$lower,
        upper = bootstrap_results[[name]]$upper,
        stringsAsFactors = FALSE
    )
}))

write.xlsx(results_df, "results_subgroup.xlsx", overwrite = TRUE)

subgroup_structure <- list(
    "Age, y" = c("<65", "≥65"),
    "Sex" = c("Male", "Female"),
    "Types of stroke" = c("Ischemic stroke", "ICH", "SAH"),
    "CKD" = c("Yes CKD", "No CKD")
)

table_text_list <- list()
plot_mean_list <- list()
plot_lower_list <- list()
plot_upper_list <- list()
is_summary_vec <- c()

table_text_list[[1]] <- c("Subgroup", "N (%)", "HR [95% CI]")
plot_mean_list[[1]] <- NA
plot_lower_list[[1]] <- NA
plot_upper_list[[1]] <- NA
is_summary_vec <- c(is_summary_vec, TRUE)

for (group_name in names(subgroup_structure)) {
    table_text_list[[length(table_text_list) + 1]] <- c(group_name, "", "")
    plot_mean_list[[length(plot_mean_list) + 1]] <- NA
    plot_lower_list[[length(plot_lower_list) + 1]] <- NA
    plot_upper_list[[length(plot_upper_list) + 1]] <- NA
    is_summary_vec <- c(is_summary_vec, TRUE)

    for (sub_name in subgroup_structure[[group_name]]) {
        if (sub_name %in% results_df$subgroup) {
            idx <- which(results_df$subgroup == sub_name)
            n_val <- n_vec[sub_name]
            p_val <- percent_vec[sub_name]

            table_text_list[[length(table_text_list) + 1]] <- c(
                paste0("  ", sub_name),
                paste0(n_val, " (", sprintf("%.1f", p_val), "%)"),
                paste0(
                    sprintf("%.2f", results_df$hr[idx]),
                    " (", sprintf("%.2f", results_df$lower[idx]),
                    "-", sprintf("%.2f", results_df$upper[idx]), ")"
                )
            )
            plot_mean_list[[length(plot_mean_list) + 1]] <- results_df$hr[idx]
            plot_lower_list[[length(plot_lower_list) + 1]] <- results_df$lower[idx]
            plot_upper_list[[length(plot_upper_list) + 1]] <- results_df$upper[idx]
            is_summary_vec <- c(is_summary_vec, FALSE)
        }
    }
}

table_text <- do.call(rbind, table_text_list)
plot_mean <- unlist(plot_mean_list)
plot_lower <- unlist(plot_lower_list)
plot_upper <- unlist(plot_upper_list)

x_min <- min(plot_lower, na.rm = TRUE)
x_max <- max(plot_upper, na.rm = TRUE)
hr_ticks <- pretty(c(x_min, x_max), n = 5)

forestplot(
    labeltext = table_text,
    mean = plot_mean,
    lower = plot_lower,
    upper = plot_upper,
    is.summary = is_summary_vec,
    graphwidth = grid::unit(0.3, "npc"),
    hrzl_lines = list("2" = grid::gpar(lwd = 1.5, col = "#444444")),
    col = fpColors(
        box = "#2b5461",
        line = "#2c2e35",
        summary = "#2b5461"
    ),
    boxsize = 0.3,
    lwd.ci = 2,
    xlog = FALSE,
    xlim = c(0.5, 1.2),
    xticks = seq(0.5, 1.2, 0.1),
    xlab = "Hazard Ratio (95% CI)",
    zero = 1,
    ci_col = "black",
    line.margin = 0.05,
    txt_gp = fpTxtGp(
        label = grid::gpar(cex = 1.2),
        ticks = grid::gpar(cex = 1.1),
        xlab = grid::gpar(cex = 1.3),
        summary = grid::gpar(fontface = "bold")
    ),
    lineheight = "auto",
    colgap = grid::unit(5, "mm"),
    align = c("l", "l", "l"),
    graph.pos = 3
)



# E-value
# If the outcome incidence < 10%, rare=TRUE; otherwise, rare=FALSE
ev_hr <- evalues.HR(
    est = 0.805884,
    lo = 0.740521,
    hi = 0.864717,
    rare = FALSE
)
ev_hr

rr_point <- ev_hr["RR", "point"]
rr_lower <- ev_hr["RR", "lower"]
rr_upper <- ev_hr["RR", "upper"]

bias_plot_overlay <- function(RR_ci, RR_pt, xmax) {
    x <- seq(0, xmax, 0.01)
    if (RR_ci < 1) RR_ci <- 1 / RR_ci
    if (RR_pt < 1) RR_pt <- 1 / RR_pt

    par(cex.lab = 1.1)
    plot(x, x,
        lty = 2, col = "white", type = "l", xaxs = "i",
        yaxs = "i", xaxt = "n", yaxt = "n",
        xlab = expression(RR[EU]),
        ylab = expression(RR[UD]), xlim = c(0, xmax),
        ylim = c(0, xmax),
        main = "5-day delirium incidence"
    )

    abline(h = seq(0, xmax, by = 1), col = "#e2e2e3", lty = 1, lwd = 1.5)

    x1 <- seq(RR_ci, xmax, by = 1e-3)
    y1 <- RR_ci * (RR_ci - 1) / (x1 - RR_ci) + RR_ci
    lines(x1, y1, type = "l", col = "#b42d34", lty = 2, lwd = 2.5)

    x2 <- seq(RR_pt, xmax, by = 1e-3)
    y2 <- RR_pt * (RR_pt - 1) / (x2 - RR_pt) + RR_pt
    lines(x2, y2, type = "l", col = "#1e4681", lwd = 2.5)

    high1 <- RR_ci + sqrt(RR_ci * (RR_ci - 1))
    points(high1, high1, pch = 19, col = "#b42d34")
    label1 <- paste("(", sprintf("%.2f", high1), ", ", sprintf("%.2f", high1), ")", sep = "")
    text(1, 1, label1, col = "#b42d34")

    high2 <- RR_pt + sqrt(RR_pt * (RR_pt - 1))
    points(high2, high2, pch = 19, col = "#1e4681")
    label2 <- paste("(", sprintf("%.2f", high2), ", ", sprintf("%.2f", high2), ")", sep = "")
    text(1, 2, label2, col = "#1e4681")

    axis(1, at = seq(0, xmax, by = 1), cex.axis = 1.0)
    axis(2, at = seq(0, xmax, by = 1), cex.axis = 1.0)

    legend("bottomleft",
        legend = c("E-value for HR (0.81)", "E-value for CI upper bound (0.86)"),
        col = c("#1e4681", "#b42d34"), lty = c(1, 2), lwd = 2.5, bty = "n"
    )
}

bias_plot_overlay(RR_ci = rr_upper, RR_pt = rr_point, xmax = 6)
