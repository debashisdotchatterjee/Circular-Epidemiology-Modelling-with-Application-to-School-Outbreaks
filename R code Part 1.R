#### ===============================================================
####  Circular Epidemiology for School Outbreaks
####  Datasets: outbreaks::influenza_england_1978_school (time series)
####            outbreaks::norovirus_derbyshire_2001_school (line list)
####  Author: (c) Debashis Chatterjee (Visva-Bharati University)
####  Last updated: Sys.time()
#### ===============================================================

#### ===============================================================
####  Circular Epidemiology: clean, robust, error-free script
####  Datasets: outbreaks::influenza_england_1978_school (time series)
####            outbreaks::norovirus_derbyshire_2001_school (line list)
####  Output: ./circular_epi_outputs/{figures, tables, sessionInfo.txt}
#### ===============================================================

## ---------- 0) Packages & Utilities ----------
ensure_pkgs <- function(pkgs) {
  to_get <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_get)) install.packages(to_get, repos = "https://cloud.r-project.org")
  invisible(lapply(pkgs, require, character.only = TRUE))
}

ensure_pkgs(c(
  "outbreaks","dplyr","tidyr","ggplot2","lubridate","glue","readr",
  "circular","CircStats","MASS","broom","purrr","stringr","forcats",
  "gt","knitr","scales","patchwork","tibble"
))

theme_set(theme_minimal(base_size = 12))

out_dir <- "circular_epi_outputs"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

save_plot <- function(p, fname, w=8, h=5, dpi=320){
  ggsave(filename = file.path(fig_dir, fname), plot = p, width = w, height = h, dpi = dpi)
}

save_base_plot <- function(fname, expr, w=800, h=800, res=150){
  png(filename = file.path(fig_dir, fname), width = w, height = h, res = res)
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(force(expr))
}

write_table <- function(df, fname_csv, title=NULL){
  readr::write_csv(df, file.path(tab_dir, fname_csv))
  if (!is.null(title)) cat("\n", title, "\n", sep = "")
  print(knitr::kable(df, digits = 4))
}

## Map Date -> weekday factor and angle on circle (P=7)
wk_levels <- c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")
weekday_factor <- function(x) {
  lubridate::wday(x, week_start = 1, label = TRUE, abbr = TRUE) |>
    as.character() |> factor(levels = wk_levels, ordered = TRUE)
}
weekday_angle <- function(x) {           # radians on [0, 2π)
  w <- as.integer(lubridate::wday(x, week_start = 1)) - 1L  # 0..6 (Mon=0)
  2*pi*w/7
}
angle_to_weekday <- function(theta){     # vectorised
  idx <- ((round((theta %% (2*pi)) / (2*pi/7))) %% 7) + 1L
  wk_levels[idx]
}

## VM helpers (robust)
vm_summary <- function(th){
  th <- th[is.finite(th)]
  if (length(th) < 2) {
    return(list(mu = NA_real_, mu_deg = NA_real_, mu_wk = NA_character_,
                Rbar = NA_real_, kappa = NA_real_))
  }
  th_circ <- circular::circular(th, units = "radians", modulo = "2pi")
  mu_hat   <- circular::mean.circular(th_circ)
  rho_hat  <- circular::rho.circular(th_circ)
  kappa_est <- suppressWarnings(tryCatch(circular::est.kappa(th_circ), error=function(e) NA))
  list(mu = as.numeric(mu_hat),
       mu_deg = as.numeric(mu_hat)*180/pi,
       mu_wk  = angle_to_weekday(as.numeric(mu_hat)),
       Rbar = as.numeric(rho_hat),
       kappa = as.numeric(kappa_est))
}

rayleigh_p <- function(th){
  th <- th[is.finite(th)]
  if (length(th) < 3) return(NA_real_)
  th_circ <- circular::circular(th, units="radians")
  out <- suppressWarnings(circular::rayleigh.test(th_circ))
  pv <- tryCatch(as.numeric(out$p.value), error=function(e) NA_real_)
  if (length(pv)==0) NA_real_ else pv
}

rao_p <- function(th){
  th <- th[is.finite(th)]
  if (length(th) < 4) return(NA_real_)
  th_circ <- circular::circular(th, units="radians")
  out <- suppressWarnings(circular::rao.spacing.test(th_circ))
  pv <- tryCatch(as.numeric(out$p.value), error=function(e) NA_real_)
  if (is.null(out$p.value) || length(pv)==0) NA_real_ else pv
}

## Safe Watson–Williams two-sample (fallback to CircStats if needed)
safe_ww_two_sample <- function(th1, th2){
  th1 <- th1[is.finite(th1)]; th2 <- th2[is.finite(th2)]
  if (length(th1) < 3 || length(th2) < 3) return(NA_real_)
  X <- c(th1, th2)
  grp <- factor(c(rep("g1", length(th1)), rep("g2", length(th2))))
  out <- tryCatch(
    circular::watson.williams.test(x = circular::circular(X, units="radians"), group = grp),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    pv <- tryCatch(as.numeric(out$p.value), error=function(e) NA_real_)
    return(if (length(pv)==0) NA_real_ else pv)
  }
  # Fallback: CircStats::watson.two.test (nonparametric)
  out2 <- tryCatch(CircStats::watson.two.test(th1, th2), error=function(e) NULL)
  if (is.null(out2)) return(NA_real_)
  pv2 <- tryCatch(as.numeric(out2$p.value), error=function(e) NA_real_)
  if (length(pv2)==0) NA_real_ else pv2
}

## ---------- 1) Load data ----------
data("influenza_england_1978_school", package = "outbreaks")
data("norovirus_derbyshire_2001_school", package = "outbreaks")

flu <- outbreaks::influenza_england_1978_school |>
  mutate(weekday = weekday_factor(date),
         theta   = weekday_angle(date))

noro <- outbreaks::norovirus_derbyshire_2001_school |>
  mutate(
    across(c(day_absent, start_illness, end_illness, day_vomiting), as.Date),
    wk_start = weekday_factor(start_illness),
    wk_vomit = weekday_factor(day_vomiting),
    th_start = weekday_angle(start_illness),
    th_vomit = weekday_angle(day_vomiting)
  )

cat("\nRows: flu=", nrow(flu), "  noro=", nrow(noro), "\n", sep="")

## ---------- 2) NOROVIRUS: Circular summaries ----------
noro_valid    <- noro |> filter(!is.na(th_start))
noro_valid_v  <- noro |> filter(!is.na(th_vomit))

sum_start <- vm_summary(noro_valid$th_start)
sum_vomit <- vm_summary(noro_valid_v$th_vomit)

tbl_noro_overall <- tibble::tibble(
  marker      = c("start_illness","day_vomiting"),
  mu_radians  = c(sum_start$mu,     sum_vomit$mu),
  mu_degrees  = c(sum_start$mu_deg, sum_vomit$mu_deg),
  mu_weekday  = c(sum_start$mu_wk,  sum_vomit$mu_wk),
  Rbar        = c(sum_start$Rbar,   sum_vomit$Rbar),
  kappa_est   = c(sum_start$kappa,  sum_vomit$kappa),
  p_Rayleigh  = c(rayleigh_p(noro_valid$th_start), rayleigh_p(noro_valid_v$th_vomit)),
  p_Rao       = c(rao_p(noro_valid$th_start),      rao_p(noro_valid_v$th_vomit))
)
write_table(tbl_noro_overall, "noro_overall_circular_summary.csv",
            "Norovirus overall circular summary (printed & saved):")

wk_counts_noro <- noro_valid |>
  count(wk_start) |>
  mutate(prop = n/sum(n))
write_table(wk_counts_noro, "noro_weekday_counts_start.csv",
            "Norovirus weekday counts (start_illness):")

## ---------- 3) NOROVIRUS: Plots ----------
save_base_plot("noro_rose_start.png", {
  th_c <- circular::circular(noro_valid$th_start, units="radians", modulo="2pi", zero=0, rotation="counter")
  circular::rose.diag(th_c, bins=7, prop=1.4, main="Norovirus: Start of illness (rose, 7 bins)")
  circular::plot.circular(th_c, stack=TRUE, cex=0.4, axes=FALSE, points.col="gray30", main="")
})

save_base_plot("noro_density_circular.png", {
  th_c <- circular::circular(noro_valid$th_start, units="radians")
  dhat <- circular::density.circular(th_c, bw = 20)
  plot(dhat, main="Norovirus: Circular kernel density (start_illness)")
})

p_noro_polar <- wk_counts_noro |>
  ggplot(aes(x = wk_start, y = n)) +
  geom_col(width = 1, alpha = 0.9) +
  coord_polar() +
  labs(title="Norovirus: Weekday distribution (start_illness)", x=NULL, y="Count") +
  theme(axis.text.y = element_blank(), panel.grid.minor = element_blank())
save_plot(p_noro_polar, "noro_polar_weekday_start.png")

save_base_plot("noro_rose_start_vs_vomit.png", {
  par(mfrow=c(1,2))
  th1 <- circular::circular(noro_valid$th_start, units="radians")
  th2 <- circular::circular(noro_valid_v$th_vomit, units="radians")
  circular::rose.diag(th1, bins=7, prop=1.4, main="Start of illness")
  circular::rose.diag(th2, bins=7, prop=1.4, main="Day of vomiting")
  par(mfrow=c(1,1))
})

ww_p <- safe_ww_two_sample(noro_valid$th_start, noro_valid_v$th_vomit)
tbl_noro_compare <- tibble::tibble(
  test = "Watson–Williams (or Watson two-sample fallback)",
  p_value = ww_p
)
write_table(tbl_noro_compare, "noro_start_vs_vomit_ww.csv",
            "Norovirus start vs vomiting: mean-direction comparison:")

## ---------- 4) INFLUENZA: Circular summaries from counts ----------
## Expand daily counts into pseudo-individual angles for circular summaries
flu_expand <- flu |>
  rowwise() |>
  mutate(angles = list(rep(theta, in_bed))) |>
  ungroup() |>
  tidyr::unnest(angles)

flu_vm <- vm_summary(flu_expand$angles)
tbl_flu_circ <- tibble::tibble(
  dataset     = "influenza_england_1978_school (in_bed)",
  mu_radians  = flu_vm$mu, mu_degrees = flu_vm$mu_deg, mu_weekday = flu_vm$mu_wk,
  Rbar        = flu_vm$Rbar, kappa_est = flu_vm$kappa,
  p_Rayleigh  = rayleigh_p(flu_expand$angles),
  p_Rao       = rao_p(flu_expand$angles)
)
write_table(tbl_flu_circ, "flu_circular_summary.csv",
            "Influenza (counts expanded) circular summary:")

wk_counts_flu <- flu |>
  group_by(weekday) |>
  summarise(count = sum(in_bed), .groups="drop") |>
  mutate(prop = count/sum(count))
write_table(wk_counts_flu, "flu_weekday_counts.csv",
            "Influenza: total in_bed by weekday:")

## ---------- 5) INFLUENZA: Harmonic GLMs ----------
flu_glm_data <- flu |>
  mutate(t_index = row_number()-1L, # t=0..13
         c1 = cos(2*pi*t_index/7), s1 = sin(2*pi*t_index/7),
         c2 = cos(4*pi*t_index/7), s2 = sin(4*pi*t_index/7))

glm0 <- glm(in_bed ~ 1, data=flu_glm_data, family=poisson())
glm1 <- glm(in_bed ~ c1 + s1, data=flu_glm_data, family=poisson())
glm2 <- glm(in_bed ~ c1 + s1 + c2 + s2, data=flu_glm_data, family=poisson())

over_phi <- function(m){
  sum(residuals(m, type="pearson")^2) / df.residual(m)
}
phis <- c(M0=over_phi(glm0), M1=over_phi(glm1), M2=over_phi(glm2))
print(phis)

use_negbin <- any(phis > 1.5)
if (use_negbin) {
  glm1_nb <- MASS::glm.nb(in_bed ~ c1 + s1, data=flu_glm_data)
  glm2_nb <- MASS::glm.nb(in_bed ~ c1 + s1 + c2 + s2, data=flu_glm_data)
  models <- list(Pois_M0=glm0, Pois_M1=glm1, Pois_M2=glm2, NB_M1=glm1_nb, NB_M2=glm2_nb)
} else {
  models <- list(Pois_M0=glm0, Pois_M1=glm1, Pois_M2=glm2)
}

model_tbl <- purrr::map_dfr(names(models), function(nm){
  m <- models[[nm]]
  tibble::tibble(
    model = nm,
    family = class(m)[1],
    k = length(coef(m)),
    AIC = AIC(m),
    BIC = BIC(m)
  )
}) |> arrange(AIC)
write_table(model_tbl, "flu_model_selection.csv", "Influenza GLM model selection (AIC/BIC):")

coef_source <- if (use_negbin) models[[if("NB_M2" %in% names(models)) "NB_M2" else "NB_M1"]] else models[["Pois_M2"]]
co <- coef(coef_source)
beta_c1 <- unname(co["c1"])
beta_s1 <- unname(co["s1"])
amp1 <- sqrt(beta_c1^2 + beta_s1^2)
mu_glm <- atan2(beta_s1, beta_c1)                     # radians
mu_glm_deg <- mu_glm*180/pi
mu_glm_wk  <- angle_to_weekday(mu_glm)

tbl_phase <- tibble::tibble(
  model_used = names(models)[which(purrr::map_lgl(models, identical, y=coef_source))],
  beta_c1 = beta_c1, beta_s1 = beta_s1,
  amplitude = amp1, phase_rad = mu_glm, phase_deg = mu_glm_deg,
  phase_weekday = mu_glm_wk
)
write_table(tbl_phase, "flu_phase_amplitude.csv",
            "Influenza fundamental harmonic: amplitude & phase:")

## ---------- 6) INFLUENZA: Plots ----------
flu_glm_best <- coef_source
flu_glm_data$fit <- predict(flu_glm_best, type="response")

p_ts <- ggplot(flu_glm_data, aes(x=date, y=in_bed)) +
  geom_line(linewidth=1) + geom_point() +
  geom_line(aes(y=fit), linetype=2) +
  labs(title="Influenza (boarding school, 1978): in_bed and harmonic fit",
       y="In bed", x=NULL, subtitle=paste0("Model: ", tbl_phase$model_used[1])) +
  theme(panel.grid.minor = element_blank())
save_plot(p_ts, "flu_time_series_fit.png")

save_base_plot("flu_rose_counts.png", {
  th_c <- circular::circular(flu_expand$angles, units="radians")
  circular::rose.diag(th_c, bins=7, prop=1.6, main="Influenza: Weekly rose (weighted by in_bed)")
})

p_flu_polar <- wk_counts_flu |>
  ggplot(aes(x = weekday, y = count)) +
  geom_col(width = 1, alpha = 0.9) +
  coord_polar() +
  labs(title="Influenza: Weekday distribution (in_bed total)",
       x=NULL, y="Total in_bed") +
  theme(axis.text.y = element_blank(), panel.grid.minor = element_blank())
save_plot(p_flu_polar, "flu_polar_weekday.png")

day_grid <- tibble::tibble(t_index = seq(0, 13, by=0.1)) |>
  mutate(c1 = cos(2*pi*t_index/7), s1 = sin(2*pi*t_index/7),
         c2 = cos(4*pi*t_index/7), s2 = sin(4*pi*t_index/7))
pred <- predict(flu_glm_best, newdata = day_grid, type="link", se.fit = TRUE)
day_grid <- day_grid |>
  mutate(link = pred$fit, se = pred$se.fit,
         mu = exp(link),
         lower = exp(link - 1.96*se),
         upper = exp(link + 1.96*se),
         phase_angle = (2*pi*(t_index %% 7))/7,
         weekday = factor(angle_to_weekday(phase_angle), levels = wk_levels, ordered = TRUE))

p_cycle <- ggplot(day_grid, aes(x = phase_angle, y = mu)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(linewidth=1) +
  scale_x_continuous(breaks = 2*pi*(0:6)/7, labels = wk_levels, limits=c(0, 2*pi)) +
  labs(title="Influenza: Fitted weekly cycle (fundamental harmonic)",
       x="Phase (Mon→Sun)", y="Mean in_bed (model scale)")
save_plot(p_cycle, "flu_weekly_cycle.png", w=9, h=4.5)

## ---------- 7) NOROVIRUS: Class-level contrasts (top 6 classes) ----------
top_classes <- noro_valid |>
  count(class, sort=TRUE) |>
  slice_head(n=6) |>
  pull(class)
noro_top <- noro_valid |> filter(class %in% top_classes)

class_summ <- noro_top |>
  group_by(class) |>
  reframe(
    n = dplyr::n(),
    mu        = vm_summary(th_start)$mu,
    mu_deg    = vm_summary(th_start)$mu_deg,
    mu_wk     = vm_summary(th_start)$mu_wk,
    Rbar      = vm_summary(th_start)$Rbar,
    kappa     = vm_summary(th_start)$kappa,
    p_Rayleigh= rayleigh_p(th_start),
    p_Rao     = rao_p(th_start)
  )
write_table(class_summ, "noro_top_classes_circular.csv",
            "Norovirus: top classes circular summaries:")

wk_class <- noro_top |>
  mutate(wk = wk_start) |>
  count(class, wk)

p_noro_class <- wk_class |>
  ggplot(aes(x = wk, y = n)) +
  geom_col(width = 1, alpha = 0.9) +
  coord_polar() +
  facet_wrap(~ class, ncol = 3) +
  labs(title="Norovirus: Weekday distribution by class (top 6)",
       x=NULL, y="Count") +
  theme(axis.text.y=element_blank(), panel.grid.minor = element_blank())
save_plot(p_noro_class, "noro_polar_by_class.png", w=10, h=7)

## ---------- 8) Save coefficient tables nicely ----------
coef_tables <- purrr::map(models, ~ broom::tidy(.x)) |>
  purrr::imap(~ dplyr::mutate(.x, model = .y)) |>
  dplyr::bind_rows() |>
  dplyr::select(model, term, estimate, std.error, statistic, p.value)
write_table(coef_tables, "flu_glm_coefficients.csv",
            "Influenza: GLM coefficients:")

## ---------- 9) Session info ----------
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))
cat("\nAll done. Outputs saved in:", normalizePath(out_dir), "\n")
