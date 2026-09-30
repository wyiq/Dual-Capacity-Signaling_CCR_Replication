# ==============================================================================
# Douyin Video State Capacity Signaling Analysis
# ==============================================================================
# This script performs the following analyses:
# Part I: Impact of signaling on Engagement
# ==============================================================================

rm(list = ls())


# Load required packages -------------------------------------------------------
library(tidyverse)
library(nnet)        # For multinomial logistic regression
library(broom)       # For tidy model outputs
library(scales)      # For formatting
library(texreg)      # For LaTeX tables
library(patchwork)   # For combining plots

# ==============================================================================
# Prompt-style plotting helpers
# ==============================================================================

if (requireNamespace("showtext", quietly = TRUE)) {
  showtext::showtext_auto()
}

theme_set(theme_minimal(base_size = 12))

theme_rr <- function() {
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
}

normalize_emm <- function(tbl) {
  nm <- names(tbl)
  
  if (all(c("lower.CL", "upper.CL") %in% nm)) {
    tbl <- tbl %>% rename(conf_low = lower.CL, conf_high = upper.CL)
  } else if (all(c("asymp.LCL", "asymp.UCL") %in% nm)) {
    tbl <- tbl %>% rename(conf_low = asymp.LCL, conf_high = asymp.UCL)
  } else if (all(c("LCL", "UCL") %in% nm)) {
    tbl <- tbl %>% rename(conf_low = LCL, conf_high = UCL)
  } else if (all(c("lower.HPD", "upper.HPD") %in% nm)) {
    tbl <- tbl %>% rename(conf_low = lower.HPD, conf_high = upper.HPD)
  }
  
  if (!all(c("conf_low", "conf_high") %in% names(tbl))) {
    stop("Could not find confidence interval columns in emmeans output.")
  }
  
  tbl
}

pick_est_col <- function(tbl, candidates) {
  out <- intersect(candidates, names(tbl))
  if (length(out) == 0) {
    stop("Could not find an estimate column. Available columns: ",
         paste(names(tbl), collapse = ", "))
  }
  out[1]
}

extract_signal_emm <- function(emm_obj, conf_level = 0.95, type = NULL) {
  emm_sum <- if (is.null(type)) {
    summary(emm_obj, infer = TRUE, level = conf_level)
  } else {
    summary(emm_obj, type = type, infer = TRUE, level = conf_level)
  }
  
  tbl <- emm_sum %>%
    as_tibble() %>%
    normalize_emm()
  
  est_col <- pick_est_col(
    tbl,
    candidates = c("response", "emmean", "estimate", "rate", "prob")
  )
  
  tbl %>%
    filter(sc_signal_3bcat %in% c("Administrative", "Coercive")) %>%
    transmute(
      signal = sc_signal_3bcat,
      pred = .data[[est_col]],
      conf_low = conf_low,
      conf_high = conf_high,
      ci_level = scales::percent(conf_level, accuracy = 1)
    )
}

extract_ratio_emm <- function(emm_obj, conf_level = 0.95) {
  lev <- emm_obj@levels$sc_signal_3bcat
  
  if (!all(c("Coercive", "Administrative") %in% lev)) {
    stop("Expected Coercive and Administrative in sc_signal_3bcat levels.")
  }
  
  contrast_vec <- rep(0, length(lev))
  names(contrast_vec) <- lev
  contrast_vec["Administrative"] <- 1
  contrast_vec["Coercive"] <- -1
  
  ratio_tbl <- contrast(
    emm_obj,
    method = list("Administrative / Coercive" = contrast_vec)
  ) %>%
    summary(type = "response", infer = TRUE, level = conf_level) %>%
    as_tibble() %>%
    normalize_emm()
  
  est_col <- pick_est_col(
    ratio_tbl,
    candidates = c("ratio", "response", "estimate", "odds.ratio")
  )
  
  ratio_tbl %>%
    transmute(
      estimate = .data[[est_col]],
      conf_low = conf_low,
      conf_high = conf_high,
      p_value = if ("p.value" %in% names(ratio_tbl)) p.value else NA_real_,
      ci_level = scales::percent(conf_level, accuracy = 1)
    )
}

format_metric_signal <- function(df) {
  df %>%
    mutate(
      metric = factor(
        metric,
        levels = c("digg", "comment", "share", "collect"),
        labels = c("Likes", "Comments", "Shares", "Collections")
      ),
      signal = factor(signal, levels = c("Coercive", "Administrative")),
      ci_level = factor(ci_level, levels = c("90%", "95%"))
    )
}

format_metric_ratio <- function(df) {
  df %>%
    mutate(
      metric = factor(
        metric,
        levels = c("digg", "comment", "share", "collect"),
        labels = c("Likes", "Comments", "Shares", "Collections")
      ),
      metric = forcats::fct_rev(metric),
      ci_level = factor(ci_level, levels = c("90%", "95%"))
    )
}

plot_pred_signal_rr <- function(plot_data,
                                output_filename,
                                title,
                                subtitle,
                                y_lab,
                                out_dir = output_dir,
                                width = 8,
                                height = 5) {
  
  plot_data <- format_metric_signal(plot_data)
  
  plot_data_90 <- plot_data %>% filter(ci_level == "90%")
  plot_data_95 <- plot_data %>% filter(ci_level == "95%")
  
  pos <- position_dodge(width = 0.45)
  
  p <- ggplot() +
    # 90% CI: thick line
    geom_linerange(
      data = plot_data_90,
      aes(
        x = signal,
        ymin = conf_low,
        ymax = conf_high,
        color = signal,
        group = signal
      ),
      linewidth = 1,
      position = pos
    ) +
    # 95% CI: thin line
    geom_linerange(
      data = plot_data_95,
      aes(
        x = signal,
        ymin = conf_low,
        ymax = conf_high,
        color = signal,
        group = signal
      ),
      linewidth = 0.5,
      position = pos
    ) +
    geom_point(
      data = plot_data_95,
      aes(x = signal, y = pred, color = signal, shape = signal),
      size = 3,
      position = pos
    ) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_color_brewer(palette = "Set1", name = "Signal") +
    scale_shape_manual(values = c(16, 17), name = "Signal") +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = y_lab,
      caption = "Thick lines: 90% CI; Thin lines: 95% CI"
    ) +
    theme_rr() +
    theme(
      axis.text.x = element_text(size = 11, hjust = 0.5)
    )
  
  ggsave(
    file.path(out_dir, output_filename),
    p,
    width = width,
    height = height,
    device = "pdf"
  )
  
  cat("Saved:", output_filename, "\n")
  return(p)
}

plot_ratio_rr <- function(plot_data,
                          output_filename,
                          title,
                          subtitle,
                          x_lab = "Ratio of mean engagement",
                          out_dir = output_dir,
                          width = 7,
                          height = 4.2,
                          color_choice = "steelblue") {
  
  plot_data <- format_metric_ratio(plot_data)
  
  plot_data_90 <- plot_data %>% filter(ci_level == "90%")
  plot_data_95 <- plot_data %>% filter(ci_level == "95%")
  
  use_log_x <- all(plot_data$estimate > 0, na.rm = TRUE) &&
    all(plot_data$conf_low > 0, na.rm = TRUE) &&
    all(plot_data$conf_high > 0, na.rm = TRUE)
  
  caption_text <- "Thick lines: 90% CI; Thin lines: 95% CI"
  if (use_log_x) {
    caption_text <- paste0(caption_text, "; horizontal axis on log scale")
  }
  
  p <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    # 90% CI: thick
    geom_linerange(
      data = plot_data_90,
      aes(y = metric, xmin = conf_low, xmax = conf_high),
      linewidth = 1.5,
      color = color_choice
    ) +
    # 95% CI: thin
    geom_linerange(
      data = plot_data_95,
      aes(y = metric, xmin = conf_low, xmax = conf_high),
      linewidth = 0.5,
      color = color_choice
    ) +
    geom_point(
      data = plot_data_95,
      aes(x = estimate, y = metric),
      size = 3,
      color = color_choice
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_lab,
      y = NULL,
      caption = caption_text
    ) +
    theme_rr() +
    theme(
      axis.text.y = element_text(size = 12)
    )
  
  if (use_log_x) {
    p <- p +
      scale_x_log10(labels = scales::label_number(accuracy = 0.01))
  }
  
  ggsave(
    file.path(out_dir, output_filename),
    p,
    width = width,
    height = height,
    device = "pdf"
  )
  
  cat("Saved:", output_filename, "\n")
  return(p)
}


# Create output directory ------------------------------------------------------
output_dir <- "Figures_RR"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load data --------------------------------------------------------------------
# Replace with your actual data file path
data <- read_rds("data/analysis.rds")

str(data)

# Add mplatform_followers_count by city_name -------------------------------

add_city_suffix <- function(x) {
  x <- stringr::str_squish(as.character(x))
  dplyr::if_else(stringr::str_ends(x, "市"), x, paste0(x, "市"))
}

followers_by_city <- tibble::tribble(
  ~city_name, ~mplatform_followers_count,
  "济南",       4360500,
  "重庆",       3295387,
  "合肥",       2299045,
  "上海",       1879717,
  "杭州",       1536406,
  "武汉",       1435782,
  "海口",       1167444,
  "天津",       1090066,
  "西安",        932176,
  "南京",        881402,
  "广州",        730837,
  "兰州",        724318,
  "长沙",        601785,
  "北京",        593762,
  "成都",        555535,
  "乌鲁木齐",     493303,
  "昆明",        460448,
  "银川",        369150,
  "福州",        322431,
  "哈尔滨",       275793,
  "沈阳",        257714,
  "南宁",        221563,
  "贵阳",        209273,
  "呼和浩特",     156906,
  "太原",        109277,
  "拉萨",        102457,
  "郑州",         69237,
  "南昌",         32897,
  "长春",         30479,
  "西宁",         27096,
  "石家庄",       16465
) %>%
  mutate(
    city_name = add_city_suffix(city_name),
    mplatform_followers_count = as.numeric(mplatform_followers_count)
  )

stopifnot(!anyDuplicated(followers_by_city$city_name))

data <- data %>%
  mutate(
    city_name = add_city_suffix(city_name)
  ) %>%
  left_join(followers_by_city, by = "city_name") %>%
  mutate(
    city_name_cn = city_name,
    l_mplatform_followers_count = log1p(mplatform_followers_count)
  )

missing_followers <- data %>%
  filter(is.na(mplatform_followers_count)) %>%
  distinct(city_name) %>%
  arrange(city_name)

if (nrow(missing_followers) > 0) {
  print(missing_followers)
  warning("Some cities did not match the follower-count table.")
} else {
  cat("mplatform_followers_count merged successfully for all observations.\n")
}

# Recode city names --------------------------------------------------------

data <- data %>%
  mutate(
    city_name = case_match(
      city_name,
      "西宁市"     ~ "Xining (西宁)",
      "海口市"     ~ "Haikou (海口)",
      "南昌市"     ~ "Nanchang (南昌)",
      "福州市"     ~ "Fuzhou (福州)",
      "贵阳市"     ~ "Guiyang (贵阳)",
      "杭州市"     ~ "Hangzhou (杭州)",
      "哈尔滨市"   ~ "Harbin (哈尔滨)",
      "乌鲁木齐市" ~ "Urumqi (乌鲁木齐)",
      "广州市"     ~ "Guangzhou (广州)",
      "昆明市"     ~ "Kunming (昆明)",
      "上海市"     ~ "Shanghai (上海)",
      "呼和浩特市" ~ "Hohhot (呼和浩特)",
      "成都市"     ~ "Chengdu (成都)",
      "西安市"     ~ "Xi'an (西安)",
      "兰州市"     ~ "Lanzhou (兰州)",
      "长沙市"     ~ "Changsha (长沙)",
      "银川市"     ~ "Yinchuan (银川)",
      "石家庄市"   ~ "Shijiazhuang (石家庄)",
      "天津市"     ~ "Tianjin (天津)",
      "南京市"     ~ "Nanjing (南京)",
      "南宁市"     ~ "Nanning (南宁)",
      "太原市"     ~ "Taiyuan (太原)",
      "重庆市"     ~ "Chongqing (重庆)",
      "合肥市"     ~ "Hefei (合肥)",
      "郑州市"     ~ "Zhengzhou (郑州)",
      "武汉市"     ~ "Wuhan (武汉)",
      "长春市"     ~ "Changchun (长春)",
      "沈阳市"     ~ "Shenyang (沈阳)",
      "济南市"     ~ "Jinan (济南)",
      "北京市"     ~ "Beijing (北京)",
      "拉萨市"     ~ "Lhasa (拉萨)",
      .default = city_name
    )
  )

unique(data$city_name)

# Make it a factor in alphabetical order
data <- data %>%
  mutate(
    city_name = factor(city_name, levels = sort(unique(city_name)))
  )

summary(data$city_name)

# Data preparation -------------------------------------------------------------
data <- data %>%
  mutate(
    # Convert categorical variables to factors
    # "No Signaling" is the reference category (first level) for intuitive interpretation
    sc_signal_4cat = factor(sc_signal_4cat, 
                            levels = c("No Signaling", "Coercive", "Administrative", 
                                       "Coercive + Administrative")),
    sc_signal_3bcat = factor(sc_signal_3bcat,
                             levels = c("No Signaling", "Coercive", "Administrative")),
    city_name = factor(city_name),
    year = factor(year),
    
    # Convert date
    create_date = as.Date(create_date)
  )

cat("Data loaded and prepared successfully.\n")
cat("Total observations:", nrow(data), "\n\n")



# ============================================
# Engagement outcomes: Admin vs Coercive (log)
# ============================================

# Packages
library(tidyverse)
library(lubridate)
library(emmeans)
library(patchwork)
theme_set(theme_minimal(base_size = 12))

# Ensure signals and FEs
data <- data %>%
  mutate(
    # Use the 3-category DV for signal; ensure expected levels/order
    sc_signal_3bcat = factor(sc_signal_3bcat,
                             levels = c("No Signaling", "Coercive", "Administrative")),
    # Year FE from create_date (if you already have a year factor you want to use, replace year_fe with it)
    year_fe = factor(lubridate::year(create_date)),
    # Log(1 + x) transforms for outcomes and duration
    l_engagement_digg     = log1p(engagement_digg),
    l_engagement_comment  = log1p(engagement_comment),
    l_engagement_share    = log1p(engagement_share),
    l_engagement_collect  = log1p(engagement_collect),
    l_duration            = log1p(duration)
  )

# Ensure signals and FEs
data <- data %>%
  mutate(
    sc_signal_3bcat = factor(
      sc_signal_3bcat,
      levels = c("No Signaling", "Coercive", "Administrative")
    ),
    year_fe = factor(lubridate::year(create_date)),
    l_engagement_digg     = log1p(engagement_digg),
    l_engagement_comment  = log1p(engagement_comment),
    l_engagement_share    = log1p(engagement_share),
    l_engagement_collect  = log1p(engagement_collect),
    l_duration            = log1p(duration)
  )

fit_and_emm <- function(df, outcome_var) {
  fml <- as.formula(
    paste0(
      outcome_var,
      " ~ sc_signal_3bcat + l_duration + has_cocreater + city_name + year_fe"
    )
  )
  
  mod <- lm(fml, data = df)
  
  emm_obj <- emmeans(
    mod,
    specs = ~ sc_signal_3bcat,
    weights = "proportional"
  )
  
  list(
    model = mod,
    emm_95 = extract_signal_emm(emm_obj, conf_level = 0.95),
    emm_90 = extract_signal_emm(emm_obj, conf_level = 0.90)
  )
}

outs <- c(
  digg    = "l_engagement_digg",
  comment = "l_engagement_comment",
  share   = "l_engagement_share",
  collect = "l_engagement_collect"
)

results <- purrr::imap(outs, ~ fit_and_emm(data, .x))

preds <- purrr::imap_dfr(
  results,
  ~ bind_rows(.x$emm_95, .x$emm_90) %>% mutate(metric = .y),
  .id = NULL
)

p_pred <- plot_pred_signal_rr(
  plot_data = preds,
  output_filename = "figure_pred_engagement_admin_vs_coercive_log.pdf",
  title = "Predicted Engagement by Signal Type",
  subtitle = "Administrative vs Coercive; controls: log(duration), has co-creator, city and year fixed effects",
  y_lab = "Predicted log(1 + engagement)",
  width = 8,
  height = 5
)

cat("Saved: figure_pred_engagement_admin_vs_coercive_log.pdf\n")

# ===========================
# Negative binomial refit (robust)
# ===========================
library(tidyverse)
library(lubridate)
library(emmeans)
library(glmmTMB)   # robust NB engine

# Prepare variables / FEs
use_date_fe <- FALSE   # set TRUE to use daily date FE (can be heavy)
engine <- "glmmTMB"    # options: "glmmTMB" (default, robust) or "glm.nb"

data <- data %>%
  mutate(
    sc_signal_3bcat = factor(sc_signal_3bcat,
                             levels = c("No Signaling", "Coercive", "Administrative")),
    year_fe = factor(lubridate::year(create_date)),
    l_duration = log1p(duration),
    city_name = factor(city_name)
  )

fe_term <- if (use_date_fe) "factor(create_date)" else "year_fe"

# Outcomes
outs <- c(
  digg    = "engagement_digg",
  comment = "engagement_comment",
  share   = "engagement_share",
  collect = "engagement_collect"
)

# Helper: unify CI and estimate column names from emmeans
normalize_emm <- function(tbl) {
  nm <- names(tbl)
  # unify CI names
  if ("lower.CL" %in% nm && "upper.CL" %in% nm) {
    tbl <- rename(tbl, conf_low = lower.CL, conf_high = upper.CL)
  } else if ("asymp.LCL" %in% nm && "asymp.UCL" %in% nm) {
    tbl <- rename(tbl, conf_low = asymp.LCL, conf_high = asymp.UCL)
  } else if ("LCL" %in% nm && "UCL" %in% nm) {
    tbl <- rename(tbl, conf_low = LCL, conf_high = UCL)
  }
  tbl
}

# ===========================
# Negative binomial refit
# ===========================

use_date_fe <- FALSE
engine <- "glmmTMB"

data <- data %>%
  mutate(
    sc_signal_3bcat = factor(
      sc_signal_3bcat,
      levels = c("No Signaling", "Coercive", "Administrative")
    ),
    year_fe = factor(lubridate::year(create_date)),
    l_duration = log1p(duration),
    city_name = factor(city_name)
  )

fe_term <- if (use_date_fe) "factor(create_date)" else "year_fe"

outs <- c(
  digg    = "engagement_digg",
  comment = "engagement_comment",
  share   = "engagement_share",
  collect = "engagement_collect"
)

fit_nb_and_emm <- function(df, outcome_var) {
  fml_txt <- paste0(
    outcome_var,
    " ~ sc_signal_3bcat + l_duration + has_cocreater + city_name + ",
    fe_term
  )
  
  fml <- as.formula(fml_txt)
  
  if (engine == "glmmTMB") {
    mod <- glmmTMB(
      fml,
      data = df,
      family = nbinom2(link = "log")
    )
  } else {
    mod <- MASS::glm.nb(
      fml,
      data = df,
      control = glm.control(maxit = 200)
    )
  }
  
  emm_obj <- emmeans(
    mod,
    specs = ~ sc_signal_3bcat,
    weights = "proportional"
  )
  
  list(
    model = mod,
    emm_95 = extract_signal_emm(emm_obj, conf_level = 0.95, type = "response"),
    emm_90 = extract_signal_emm(emm_obj, conf_level = 0.90, type = "response"),
    ratio_95 = extract_ratio_emm(emm_obj, conf_level = 0.95),
    ratio_90 = extract_ratio_emm(emm_obj, conf_level = 0.90)
  )
}

results_nb <- purrr::imap(outs, ~ fit_nb_and_emm(data, .x))

preds_nb <- purrr::imap_dfr(
  results_nb,
  ~ bind_rows(.x$emm_95, .x$emm_90) %>% mutate(metric = .y),
  .id = NULL
)

p_nb_pred <- plot_pred_signal_rr(
  plot_data = preds_nb,
  output_filename = "figure_nb_pred_engagement_admin_vs_coercive.pdf",
  title = "Predicted Engagement by Signal Type",
  subtitle = "Negative binomial models; Administrative vs Coercive Signaling",
  y_lab = "Predicted mean engagement count",
  width = 8,
  height = 5
)

ratio_nb <- purrr::imap_dfr(
  results_nb,
  ~ bind_rows(.x$ratio_95, .x$ratio_90) %>% mutate(metric = .y),
  .id = NULL
)

p_nb_ratio <- plot_ratio_rr(
  plot_data = ratio_nb,
  output_filename = "figure_nb_ratio_admin_over_coercive.pdf",
  title = "Administrative / Coercive Ratio of Predicted Engagement",
  subtitle = "Negative binomial models; averaged over observed city/year mix",
  x_lab = "Administrative / Coercive ratio",
  width = 7,
  height = 4.2,
  color_choice = "steelblue"
)

# ===========================================
# LaTeX table: NB models (Admin vs Coercive)
# ===========================================
library(texreg)

cat("\n--- Creating LaTeX table for NB engagement models ---\n")

# 1. Collect the four NB models (Diggs, Comments, Shares, Collections)
nb_models <- lapply(results_nb, `[[`, "model")
names(nb_models) <- c("Diggs", "Comments", "Shares", "Collections")

# 2. Convert each glmmTMB model into a texreg object
nb_texreg_list <- lapply(nb_models, function(m) {
  s <- summary(m)
  
  # Conditional (mean) part coefficients
  coef_mat <- s$coefficients$cond
  coef_vals <- coef_mat[, "Estimate"]
  se_vals   <- coef_mat[, "Std. Error"]
  p_vals    <- coef_mat[, "Pr(>|z|)"]
  coef_names <- rownames(coef_mat)
  
  # Goodness-of-fit statistics
  gof_vals     <- c(AIC(m), as.numeric(logLik(m)), nobs(m))
  gof_names    <- c("AIC", "Log Likelihood", "Num. obs.")
  gof_decimals <- c(TRUE, TRUE, FALSE)
  
  createTexreg(
    coef.names  = coef_names,
    coef        = coef_vals,
    se          = se_vals,
    pvalues     = p_vals,
    gof.names   = gof_names,
    gof         = gof_vals,
    gof.decimal = gof_decimals
  )
})

# 3. Custom coefficient labels (shared across models)
#    Get raw names from the first model
all_coef_names_nb <- nb_texreg_list[[1]]@coef.names

custom_coef_names_nb <- all_coef_names_nb
names(custom_coef_names_nb) <- all_coef_names_nb

# Rename main coefficients of interest
if ("(Intercept)" %in% all_coef_names_nb) {
  custom_coef_names_nb["(Intercept)"] <- "(Intercept)"
}
if ("sc_signal_3bcatCoercive" %in% all_coef_names_nb) {
  custom_coef_names_nb["sc_signal_3bcatCoercive"] <-
    "Signal: Coercive (vs No Signaling)"
}
if ("sc_signal_3bcatAdministrative" %in% all_coef_names_nb) {
  custom_coef_names_nb["sc_signal_3bcatAdministrative"] <-
    "Signal: Administrative (vs No Signaling)"
}
if ("l_duration" %in% all_coef_names_nb) {
  custom_coef_names_nb["l_duration"] <- "log(1 + Duration)"
}
if ("has_cocreater" %in% all_coef_names_nb) {
  custom_coef_names_nb["has_cocreater"] <- "Has co-creator"
}

# 4. Create LaTeX table
texreg(
  nb_texreg_list,
  file = file.path(output_dir, "table_nb_engagement_admin_vs_coercive.tex"),
  custom.model.names = names(nb_models),
  custom.coef.names  = custom_coef_names_nb,
  # Drop city and year FE coefficients from the printed table
  omit.coef = "city_name|year_fe",
  caption = "Negative binomial models of Douyin engagement by signal type (Administrative vs Coercive)",
  label   = "tab:nb_engagement_admin_vs_coercive",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = "Standard errors in parentheses. %stars. Negative binomial with log link. Controls: log(1 + duration), has co-creator, city and year fixed effects.",
  booktabs = TRUE,
  dcolumn  = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "City FE" = rep("Yes", length(nb_texreg_list)),
    "Year FE" = rep("Yes", length(nb_texreg_list))
  )
)

cat("Saved: table_nb_engagement_admin_vs_coercive.tex\n")


# ==============================================================================
# Additional models with mplatform_followers_count
# ==============================================================================

cat("\n--- Refitting models with mplatform_followers_count ---\n")

followers_model_note <- c(
  "NOTE for *_with_followers outputs:",
  "Account_followers_count is merged at the city/account level.",
  "The models use log(1 + Account_followers_count).",
  "Because follower count is constant within city_name, it is perfectly collinear with city fixed effects.",
  "Therefore, the follower-adjusted models below retain year FE but omit city FE.",
  "If city FE are kept, Account_followers_count is not separately identified."
)

cat(paste(followers_model_note, collapse = "\n"), "\n")

writeLines(
  followers_model_note,
  con = file.path(output_dir, "note_models_with_followers.txt")
)

data_followers <- data %>%
  mutate(
    sc_signal_3bcat = factor(
      sc_signal_3bcat,
      levels = c("No Signaling", "Coercive", "Administrative")
    ),
    create_date = as.Date(create_date),
    year_fe = factor(lubridate::year(create_date)),
    l_engagement_digg    = log1p(engagement_digg),
    l_engagement_comment = log1p(engagement_comment),
    l_engagement_share   = log1p(engagement_share),
    l_engagement_collect = log1p(engagement_collect),
    l_duration = log1p(duration),
    mplatform_followers_count = as.numeric(mplatform_followers_count),
    l_mplatform_followers_count = log1p(mplatform_followers_count)
  ) %>%
  filter(!is.na(l_mplatform_followers_count))

normalize_emm <- function(tbl) {
  nm <- names(tbl)
  
  if ("lower.CL" %in% nm && "upper.CL" %in% nm) {
    tbl <- rename(tbl, conf_low = lower.CL, conf_high = upper.CL)
  } else if ("asymp.LCL" %in% nm && "asymp.UCL" %in% nm) {
    tbl <- rename(tbl, conf_low = asymp.LCL, conf_high = asymp.UCL)
  } else if ("LCL" %in% nm && "UCL" %in% nm) {
    tbl <- rename(tbl, conf_low = LCL, conf_high = UCL)
  }
  
  tbl
}

make_coef_map_followers <- function(coef_names) {
  keep <- !stringr::str_detect(
    coef_names,
    "^year_fe|^factor\\(create_date\\)|^city_name"
  )
  
  raw <- coef_names[keep]
  lab <- raw
  
  lab[raw == "(Intercept)"] <- "(Intercept)"
  lab[raw == "sc_signal_3bcatCoercive"] <-
    "Signal: Coercive (vs No Signaling)"
  lab[raw == "sc_signal_3bcatAdministrative"] <-
    "Signal: Administrative (vs No Signaling)"
  lab[raw == "l_duration"] <- "log(1 + Duration)"
  lab[raw %in% c("has_cocreater", "has_cocreaterTRUE", "has_cocreater1")] <-
    "Has co-creator"
  lab[raw == "l_Account_followers_count"] <-
    "log(1 + Platform followers)"
  
  stats::setNames(as.list(lab), raw)
}

# ============================================
# Log OLS engagement models with followers
# ============================================

fit_and_emm_followers <- function(df, outcome_var) {
  fml <- as.formula(
    paste0(
      outcome_var,
      " ~ sc_signal_3bcat + l_duration + has_cocreater + ",
      "l_mplatform_followers_count + year_fe"
    )
  )
  
  mod <- lm(fml, data = df)
  
  emm_obj <- emmeans(
    mod,
    specs = ~ sc_signal_3bcat,
    weights = "proportional"
  )
  
  list(
    model = mod,
    emm_95 = extract_signal_emm(emm_obj, conf_level = 0.95),
    emm_90 = extract_signal_emm(emm_obj, conf_level = 0.90)
  )
}

outs_log_followers <- c(
  digg    = "l_engagement_digg",
  comment = "l_engagement_comment",
  share   = "l_engagement_share",
  collect = "l_engagement_collect"
)

results_log_followers <- purrr::imap(
  outs_log_followers,
  ~ fit_and_emm_followers(data_followers, .x)
)

preds_log_followers <- purrr::imap_dfr(
  results_log_followers,
  ~ bind_rows(.x$emm_95, .x$emm_90) %>% mutate(metric = .y),
  .id = NULL
)

p_pred_followers <- plot_pred_signal_rr(
  plot_data = preds_log_followers,
  output_filename = "figure_pred_engagement_admin_vs_coercive_log_with_followers.pdf",
  title = "Predicted Engagement with Platform Followers",
  subtitle = "Administrative vs Coercive; follower count control and year fixed effects; city FE omitted",
  y_lab = "Predicted log(1 + engagement)",
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "figure_pred_engagement_admin_vs_coercive_log_with_followers.pdf"),
  p_pred_followers,
  width = 8,
  height = 5,
  device = "pdf"
)

cat("Saved: figure_pred_engagement_admin_vs_coercive_log_with_followers.pdf\n")

lm_models_followers <- lapply(results_log_followers, `[[`, "model")
names(lm_models_followers) <- c("Likes", "Comments", "Shares", "Collections")

texreg(
  lm_models_followers,
  file = file.path(output_dir, "table_lm_engagement_admin_vs_coercive_with_followers.tex"),
  custom.model.names = names(lm_models_followers),
  custom.coef.map = make_coef_map_followers(names(coef(lm_models_followers[[1]]))),
  caption = "OLS models of logged Douyin engagement with platform follower count",
  label = "tab:lm_engagement_admin_vs_coercive_with_followers",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = paste(
    "Standard errors in parentheses. %stars.",
    "Follower count enters as log(1 + Account_followers_count).",
    "City fixed effects are omitted because follower count is constant within city and therefore collinear with city FE.",
    "Year fixed effects are included."
  ),
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "Platform followers control" = rep("Yes", length(lm_models_followers)),
    "City FE" = rep("No", length(lm_models_followers)),
    "Year FE" = rep("Yes", length(lm_models_followers))
  )
)

cat("Saved: table_lm_engagement_admin_vs_coercive_with_followers.tex\n")

# ============================================
# Negative binomial models with followers
# ============================================

engine_followers <- "glmmTMB"

fe_term_followers <- "year_fe"

fit_nb_and_emm_followers <- function(df, outcome_var) {
  fml_txt <- paste0(
    outcome_var,
    " ~ sc_signal_3bcat + l_duration + has_cocreater + ",
    "l_mplatform_followers_count + ",
    fe_term_followers
  )
  
  fml <- as.formula(fml_txt)
  
  if (engine_followers == "glmmTMB") {
    mod <- glmmTMB(
      fml,
      data = df,
      family = nbinom2(link = "log")
    )
  } else {
    mod <- MASS::glm.nb(
      fml,
      data = df,
      control = glm.control(maxit = 200)
    )
  }
  
  emm_obj <- emmeans(
    mod,
    specs = ~ sc_signal_3bcat,
    weights = "proportional"
  )
  
  list(
    model = mod,
    emm_95 = extract_signal_emm(emm_obj, conf_level = 0.95, type = "response"),
    emm_90 = extract_signal_emm(emm_obj, conf_level = 0.90, type = "response"),
    ratio_95 = extract_ratio_emm(emm_obj, conf_level = 0.95),
    ratio_90 = extract_ratio_emm(emm_obj, conf_level = 0.90)
  )
}

outs_count_followers <- c(
  digg    = "engagement_digg",
  comment = "engagement_comment",
  share   = "engagement_share",
  collect = "engagement_collect"
)

results_nb_followers <- purrr::imap(
  outs_count_followers,
  ~ fit_nb_and_emm_followers(data_followers, .x)
)

preds_nb_followers <- purrr::imap_dfr(
  results_nb_followers,
  ~ bind_rows(.x$emm_95, .x$emm_90) %>% mutate(metric = .y),
  .id = NULL
)

p_nb_pred_followers <- plot_pred_signal_rr(
  plot_data = preds_nb_followers,
  output_filename = "figure_nb_pred_engagement_admin_vs_coercive_with_followers.pdf",
  title = "Predicted Engagement with Platform Followers",
  subtitle = "Negative binomial models; follower count control and year fixed effects; city FE omitted",
  y_lab = "Predicted mean engagement count",
  width = 8,
  height = 5
)

ratio_nb_followers <- purrr::imap_dfr(
  results_nb_followers,
  ~ bind_rows(.x$ratio_95, .x$ratio_90) %>% mutate(metric = .y),
  .id = NULL
)

p_nb_ratio_followers <- plot_ratio_rr(
  plot_data = ratio_nb_followers,
  output_filename = "figure_nb_ratio_admin_over_coercive_with_followers.pdf",
  title = "Administrative / Coercive Ratio with Platform Followers",
  subtitle = "Negative binomial models; follower count control and year fixed effects; city FE omitted",
  x_lab = "Administrative / Coercive ratio",
  width = 7,
  height = 4.2,
  color_choice = "steelblue"
)

cat("Saved: figure_nb_ratio_admin_over_coercive_with_followers.pdf\n")


# ===========================================
# LaTeX table: NB models with followers
# ===========================================

cat("\n--- Creating LaTeX table for NB models with followers ---\n")

nb_models_followers <- lapply(results_nb_followers, `[[`, "model")
names(nb_models_followers) <- c("Diggs", "Comments", "Shares", "Collections")

nb_texreg_list_followers <- lapply(nb_models_followers, function(m) {
  s <- summary(m)
  
  coef_mat <- s$coefficients$cond
  coef_vals <- coef_mat[, "Estimate"]
  se_vals   <- coef_mat[, "Std. Error"]
  p_vals    <- coef_mat[, "Pr(>|z|)"]
  coef_names <- rownames(coef_mat)
  
  gof_vals     <- c(AIC(m), as.numeric(logLik(m)), nobs(m))
  gof_names    <- c("AIC", "Log Likelihood", "Num. obs.")
  gof_decimals <- c(TRUE, TRUE, FALSE)
  
  createTexreg(
    coef.names  = coef_names,
    coef        = coef_vals,
    se          = se_vals,
    pvalues     = p_vals,
    gof.names   = gof_names,
    gof         = gof_vals,
    gof.decimal = gof_decimals
  )
})

texreg(
  nb_texreg_list_followers,
  file = file.path(output_dir, "table_nb_engagement_admin_vs_coercive_with_followers.tex"),
  custom.model.names = names(nb_models_followers),
  custom.coef.map = make_coef_map_followers(nb_texreg_list_followers[[1]]@coef.names),
  caption = "Negative binomial models of Douyin engagement with platform follower count",
  label = "tab:nb_engagement_admin_vs_coercive_with_followers",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = paste(
    "Standard errors in parentheses. %stars.",
    "Negative binomial models with log link.",
    "Follower count enters as log(1 + Account_followers_count).",
    "City fixed effects are omitted because follower count is constant within city and therefore collinear with city FE.",
    "Year fixed effects are included."
  ),
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "Platform followers control" = rep("Yes", length(nb_texreg_list_followers)),
    "City FE" = rep("No", length(nb_texreg_list_followers)),
    "Year FE" = rep("Yes", length(nb_texreg_list_followers))
  )
)

cat("Saved: table_nb_engagement_admin_vs_coercive_with_followers.tex\n")

# ==============================================================================
# Additional models: Platform followers + City FE
# ==============================================================================

cat("\n--- Refitting models with platform followers AND city FE ---\n")

data_followers_cityfe <- data %>%
  mutate(
    sc_signal_3bcat = factor(
      sc_signal_3bcat,
      levels = c("No Signaling", "Coercive", "Administrative")
    ),
    create_date = as.Date(create_date),
    year_fe = factor(lubridate::year(create_date)),
    city_name = factor(city_name),
    
    l_engagement_digg    = log1p(engagement_digg),
    l_engagement_comment = log1p(engagement_comment),
    l_engagement_share   = log1p(engagement_share),
    l_engagement_collect = log1p(engagement_collect),
    
    l_duration = log1p(duration),
    mplatform_followers_count = as.numeric(mplatform_followers_count),
    l_mplatform_followers_count = log1p(mplatform_followers_count)
  ) %>%
  filter(!is.na(l_mplatform_followers_count)) %>%
  droplevels()

# Check whether followers vary within city
city_follower_check <- data_followers_cityfe %>%
  group_by(city_name) %>%
  summarise(
    n_follower_values = n_distinct(l_mplatform_followers_count),
    .groups = "drop"
  )

if (all(city_follower_check$n_follower_values == 1)) {
  warning(
    "Account_followers_count is constant within city_name. ",
    "With city fixed effects, the follower coefficient is not separately identified. ",
    "City FE absorb this city-level control."
  )
}

# Put follower variable after city FE so that if a redundant term is dropped,
# the follower variable is more likely to be absorbed/dropped.
rhs_followers_cityfe <- c(
  "sc_signal_3bcat",
  "l_duration",
  "has_cocreater",
  "city_name",
  "year_fe",
  "l_mplatform_followers_count"
)

fit_lm_and_emm_followers_cityfe <- function(df, outcome_var) {
  
  fml <- reformulate(
    rhs_followers_cityfe,
    response = outcome_var
  )
  
  mod <- lm(fml, data = df)
  
  dropped <- names(coef(mod))[is.na(coef(mod))]
  if (length(dropped) > 0) {
    message(
      outcome_var,
      " - dropped because of collinearity: ",
      paste(dropped, collapse = ", ")
    )
  }
  
  emm_obj <- emmeans(
    mod,
    specs = ~ sc_signal_3bcat,
    weights = "proportional"
  )
  
  list(
    model = mod,
    emm_95 = extract_signal_emm(emm_obj, conf_level = 0.95),
    emm_90 = extract_signal_emm(emm_obj, conf_level = 0.90)
  )
}

outs_log_followers_cityfe <- c(
  digg    = "l_engagement_digg",
  comment = "l_engagement_comment",
  share   = "l_engagement_share",
  collect = "l_engagement_collect"
)

results_log_followers_cityfe <- purrr::imap(
  outs_log_followers_cityfe,
  ~ fit_lm_and_emm_followers_cityfe(data_followers_cityfe, .x)
)

preds_log_followers_cityfe <- purrr::imap_dfr(
  results_log_followers_cityfe,
  ~ bind_rows(.x$emm_95, .x$emm_90) %>% mutate(metric = .y),
  .id = NULL
)

p_pred_followers_cityfe <- plot_pred_signal_rr(
  plot_data = preds_log_followers_cityfe,
  output_filename = "figure_pred_engagement_admin_vs_coercive_log_with_followers_cityfe.pdf",
  title = "Predicted Engagement with Platform Followers and City FE",
  subtitle = "Administrative vs Coercive; controls include city FE, year FE, duration, co-creator, and platform followers",
  y_lab = "Predicted log(1 + engagement)",
  width = 8,
  height = 5
)

fit_nb_and_emm_followers_cityfe <- function(df, outcome_var) {
  
  fml <- reformulate(
    rhs_followers_cityfe,
    response = outcome_var
  )
  
  mod <- glmmTMB(
    fml,
    data = df,
    family = nbinom2(link = "log"),
    control = glmmTMBControl(rank_check = "adjust")
  )
  
  emm_obj <- emmeans(
    mod,
    specs = ~ sc_signal_3bcat,
    weights = "proportional"
  )
  
  list(
    model = mod,
    emm_95 = extract_signal_emm(emm_obj, conf_level = 0.95, type = "response"),
    emm_90 = extract_signal_emm(emm_obj, conf_level = 0.90, type = "response"),
    ratio_95 = extract_ratio_emm(emm_obj, conf_level = 0.95),
    ratio_90 = extract_ratio_emm(emm_obj, conf_level = 0.90)
  )
}

outs_count_followers_cityfe <- c(
  digg    = "engagement_digg",
  comment = "engagement_comment",
  share   = "engagement_share",
  collect = "engagement_collect"
)

results_nb_followers_cityfe <- purrr::imap(
  outs_count_followers_cityfe,
  ~ fit_nb_and_emm_followers_cityfe(data_followers_cityfe, .x)
)

preds_nb_followers_cityfe <- purrr::imap_dfr(
  results_nb_followers_cityfe,
  ~ bind_rows(.x$emm_95, .x$emm_90) %>% mutate(metric = .y),
  .id = NULL
)

p_nb_pred_followers_cityfe <- plot_pred_signal_rr(
  plot_data = preds_nb_followers_cityfe,
  output_filename = "figure_nb_pred_engagement_admin_vs_coercive_with_followers_cityfe.pdf",
  title = "Predicted Engagement with Platform Followers and City FE",
  subtitle = "Negative binomial models; controls include city FE, year FE, duration, co-creator, and platform followers",
  y_lab = "Predicted mean engagement count",
  width = 8,
  height = 5
)

ratio_nb_followers_cityfe <- purrr::imap_dfr(
  results_nb_followers_cityfe,
  ~ bind_rows(.x$ratio_95, .x$ratio_90) %>% mutate(metric = .y),
  .id = NULL
)

p_nb_ratio_followers_cityfe <- plot_ratio_rr(
  plot_data = ratio_nb_followers_cityfe,
  output_filename = "figure_nb_ratio_admin_over_coercive_with_followers_cityfe.pdf",
  title = "Administrative / Coercive Ratio with Platform Followers and City FE",
  subtitle = "Negative binomial models; city FE absorb city-level follower differences",
  x_lab = "Administrative / Coercive ratio",
  width = 7,
  height = 4.2,
  color_choice = "steelblue"
)

custom.note = paste(
  "Standard errors in parentheses. %stars.",
  "Models include year fixed effects and city fixed effects.",
  "Platform follower count is city-level and therefore collinear with city fixed effects;",
  "the follower control is absorbed by city FE and is not separately identified."
)

custom.gof.rows = list(
  "Platform followers control" = rep("Absorbed by City FE", 4),
  "City FE" = rep("Yes", 4),
  "Year FE" = rep("Yes", 4)
)


# ==============================================================================
# Tables: Negative binomial models of Douyin engagement by 4-category signaling
#   Table 1: without account follower control
#   Table 2: with account follower control
# ==============================================================================

library(tidyverse)
library(lubridate)
library(glmmTMB)
library(texreg)

if (!exists("output_dir")) output_dir <- "Figures_RR"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Prepare variables
# ------------------------------------------------------------------------------

# Detect follower-count variable from your data
follower_candidates <- c(
  "account_followers_count",
  "Account_followers_count",
  "mplatform_followers_count"
)

follower_var <- intersect(follower_candidates, names(data))

if (length(follower_var) == 0) {
  stop(
    "No follower-count variable found. Expected one of: ",
    paste(follower_candidates, collapse = ", ")
  )
}

follower_var <- follower_var[1]
message("Using follower-count variable: ", follower_var)

data_nb_4cat <- data %>%
  mutate(
    sc_signal_4cat = factor(
      sc_signal_4cat,
      levels = c(
        "No Signaling",
        "Coercive",
        "Administrative",
        "Coercive + Administrative"
      )
    ),
    create_date = as.Date(create_date),
    year_fe = factor(lubridate::year(create_date)),
    city_name = factor(city_name),
    l_duration = log1p(duration),
    account_followers_count = as.numeric(.data[[follower_var]]),
    l_account_followers_count = log1p(account_followers_count)
  ) %>%
  droplevels()

# Engagement outcomes
engagement_outcomes <- c(
  "Likes"       = "engagement_digg",
  "Comments"   = "engagement_comment",
  "Shares"     = "engagement_share",
  "Collections" = "engagement_collect"
)

# Follower-complete sample
data_nb_4cat_followers <- data_nb_4cat %>%
  filter(!is.na(l_account_followers_count)) %>%
  droplevels()

if (nrow(data_nb_4cat_followers) == 0) {
  stop("No observations have non-missing follower counts.")
}

# ------------------------------------------------------------------------------
# 2. Decide whether city FE can be included with follower control
# ------------------------------------------------------------------------------

city_follower_check <- data_nb_4cat_followers %>%
  group_by(city_name) %>%
  summarise(
    n_follower_values = n_distinct(l_account_followers_count, na.rm = TRUE),
    .groups = "drop"
  )

include_city_fe_with_followers <- any(city_follower_check$n_follower_values > 1)

if (!include_city_fe_with_followers) {
  message(
    "Follower count is constant within city_name. ",
    "Follower-control models will omit city FE because the follower control ",
    "is collinear with city fixed effects."
  )
}

# If you want to force city FE in follower models, uncomment:
# include_city_fe_with_followers <- TRUE

# ------------------------------------------------------------------------------
# 3. Model-fitting helper
# ------------------------------------------------------------------------------

fit_nb_4cat <- function(df,
                        outcome_var,
                        include_followers = FALSE,
                        include_city_fe = TRUE) {
  
  rhs <- c(
    "sc_signal_4cat",
    "l_duration",
    "has_cocreater"
  )
  
  if (include_followers) {
    rhs <- c(rhs, "l_account_followers_count")
  }
  
  if (include_city_fe) {
    rhs <- c(rhs, "city_name")
  }
  
  rhs <- c(rhs, "year_fe")
  
  fml <- reformulate(rhs, response = outcome_var)
  
  glmmTMB(
    fml,
    data = df,
    family = nbinom2(link = "log"),
    control = glmmTMBControl(rank_check = "adjust")
  )
}

# ------------------------------------------------------------------------------
# 4. Fit models: without follower control
# ------------------------------------------------------------------------------

nb_4cat_no_followers <- purrr::map(
  engagement_outcomes,
  ~ fit_nb_4cat(
    df = data_nb_4cat,
    outcome_var = .x,
    include_followers = FALSE,
    include_city_fe = TRUE
  )
)

names(nb_4cat_no_followers) <- names(engagement_outcomes)

# ------------------------------------------------------------------------------
# 5. Fit models: with follower control
# ------------------------------------------------------------------------------

nb_4cat_with_followers <- purrr::map(
  engagement_outcomes,
  ~ fit_nb_4cat(
    df = data_nb_4cat_followers,
    outcome_var = .x,
    include_followers = TRUE,
    include_city_fe = include_city_fe_with_followers
  )
)

names(nb_4cat_with_followers) <- names(engagement_outcomes)

# ------------------------------------------------------------------------------
# 6. texreg extraction helper for glmmTMB models
# ------------------------------------------------------------------------------

glmmTMB_to_texreg <- function(m) {
  
  s <- summary(m)
  
  coef_mat <- s$coefficients$cond
  
  coef_vals <- coef_mat[, "Estimate"]
  se_vals   <- coef_mat[, "Std. Error"]
  p_vals    <- coef_mat[, "Pr(>|z|)"]
  coef_names <- rownames(coef_mat)
  
  gof_vals <- c(
    AIC(m),
    as.numeric(logLik(m)),
    nobs(m)
  )
  
  gof_names <- c(
    "AIC",
    "Log Likelihood",
    "Num. obs."
  )
  
  gof_decimals <- c(TRUE, TRUE, FALSE)
  
  createTexreg(
    coef.names  = coef_names,
    coef        = coef_vals,
    se          = se_vals,
    pvalues     = p_vals,
    gof.names   = gof_names,
    gof         = gof_vals,
    gof.decimal = gof_decimals
  )
}

# ------------------------------------------------------------------------------
# 7. Coefficient labels for 4-category signaling
# ------------------------------------------------------------------------------

make_coef_map_4cat <- function(coef_names) {
  
  keep <- !stringr::str_detect(
    coef_names,
    "^city_name|^year_fe|^factor\\(create_date\\)"
  )
  
  raw <- coef_names[keep]
  lab <- raw
  
  lab[raw == "(Intercept)"] <- "(Intercept)"
  
  lab[raw == "sc_signal_4catCoercive"] <-
    "Signal: Coercive (vs No Signaling)"
  
  lab[raw == "sc_signal_4catAdministrative"] <-
    "Signal: Administrative (vs No Signaling)"
  
  lab[raw == "sc_signal_4catCoercive + Administrative"] <-
    "Signal: Coercive + Administrative (vs No Signaling)"
  
  lab[raw == "l_duration"] <- "log(1 + Duration)"
  
  lab[stringr::str_detect(raw, "^has_cocreater")] <-
    "Has co-creator"
  
  lab[raw == "l_account_followers_count"] <-
    "log(1 + Account followers)"
  
  stats::setNames(as.list(lab), raw)
}

# ------------------------------------------------------------------------------
# 8. Convert models to texreg objects
# ------------------------------------------------------------------------------

nb_4cat_tr_no_followers <- purrr::map(
  nb_4cat_no_followers,
  glmmTMB_to_texreg
)

nb_4cat_tr_with_followers <- purrr::map(
  nb_4cat_with_followers,
  glmmTMB_to_texreg
)

# ------------------------------------------------------------------------------
# 9. Table 1: NB models without account follower control
# ------------------------------------------------------------------------------

texreg(
  nb_4cat_tr_no_followers,
  file = file.path(output_dir, "table_nb_engagement_4cat_no_followers.tex"),
  custom.model.names = names(nb_4cat_no_followers),
  custom.coef.map = make_coef_map_4cat(
    nb_4cat_tr_no_followers[[1]]@coef.names
  ),
  caption = paste(
    "Negative binomial models of Douyin engagement by",
    "four-category state capacity signaling"
  ),
  label = "tab:nb_engagement_4cat_no_followers",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = paste(
    "Standard errors in parentheses. %stars.",
    "Coefficients are negative-binomial log-count coefficients;",
    "exponentiating a coefficient gives the incidence-rate ratio.",
    "Models include log(1 + duration), has co-creator, city fixed effects,",
    "and year fixed effects."
  ),
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "Account follower control" = rep("No", length(nb_4cat_tr_no_followers)),
    "City FE" = rep("Yes", length(nb_4cat_tr_no_followers)),
    "Year FE" = rep("Yes", length(nb_4cat_tr_no_followers))
  )
)

cat(
  "Saved:",
  file.path(output_dir, "table_nb_engagement_4cat_no_followers.tex"),
  "\n"
)

# ------------------------------------------------------------------------------
# 10. Table 2: NB models with account follower control
# ------------------------------------------------------------------------------

city_fe_label_followers <- ifelse(
  include_city_fe_with_followers,
  "Yes",
  "No"
)

city_fe_note_followers <- ifelse(
  include_city_fe_with_followers,
  "City fixed effects are included.",
  paste(
    "City fixed effects are omitted because account follower count is constant",
    "within city and is therefore collinear with city fixed effects."
  )
)

texreg(
  nb_4cat_tr_with_followers,
  file = file.path(output_dir, "table_nb_engagement_4cat_with_followers.tex"),
  custom.model.names = names(nb_4cat_with_followers),
  custom.coef.map = make_coef_map_4cat(
    nb_4cat_tr_with_followers[[1]]@coef.names
  ),
  caption = paste(
    "Negative binomial models of Douyin engagement by",
    "four-category state capacity signaling with account follower control"
  ),
  label = "tab:nb_engagement_4cat_with_followers",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = paste(
    "Standard errors in parentheses. %stars.",
    "Coefficients are negative-binomial log-count coefficients;",
    "exponentiating a coefficient gives the incidence-rate ratio.",
    "Models include log(1 + duration), has co-creator,",
    "log(1 + account followers), and year fixed effects.",
    city_fe_note_followers
  ),
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "Account follower control" = rep("Yes", length(nb_4cat_tr_with_followers)),
    "City FE" = rep(city_fe_label_followers, length(nb_4cat_tr_with_followers)),
    "Year FE" = rep("Yes", length(nb_4cat_tr_with_followers))
  )
)

cat(
  "Saved:",
  file.path(output_dir, "table_nb_engagement_4cat_with_followers.tex"),
  "\n"
)
