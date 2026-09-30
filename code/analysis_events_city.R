# ==============================================================================
# Douyin Video State Capacity Signaling Analysis
# Complete Analysis with Cluster-Robust Standard Errors
# ==============================================================================
# This script performs three main analyses:
# Part I: Prevalence of different types of state capacity signaling
# Part II: Impact of events (temporal shocks) on signaling
#          - Analysis with both naive and cluster-robust SEs (clustered by city)
#          - Both 3-category and 4-category outcomes
# Part III: Impact of city-level indicators on signaling
#          - Analysis with both naive and cluster-robust SEs (clustered by city-year)
#          - Both 3-category and 4-category outcomes
# ==============================================================================

rm(list = ls())

# Load required packages -------------------------------------------------------
library(tidyverse)
library(nnet)        # For multinomial logistic regression
library(broom)       # For tidy model outputs
library(scales)      # For formatting
library(texreg)      # For LaTeX tables
library(patchwork)   # For combining plots
library(showtext)
showtext::showtext_auto()

# Set theme for all plots
theme_set(theme_minimal(base_size = 12))

# Create output directory ------------------------------------------------------
output_dir <- "Figures_RR"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load data --------------------------------------------------------------------
data <- read_rds("data/analysis.rds")

cat("Data loaded successfully.\n")
cat("Total observations:", nrow(data), "\n\n")


data <- data %>%
  # ChineseNewYear --> Chinese_New_Year
  rename_with(~ str_replace(., "ChineseNewYear", "Chinese_New_Year"),
              starts_with("t_")) %>%
  # TwoSessions --> Two_Sessions
  rename_with(~ str_replace(., "TwoSessions", "Two_Sessions"),
              starts_with("t_")) %>%
  # NationalDay --> National_Day
  rename_with(~ str_replace(., "NationalDay", "National_Day"),
              starts_with("t_")) %>%
  # PartyCongress --> Party_Congress
  rename_with(~ str_replace(., "PartyCongress", "Party_Congress"),
              starts_with("t_")) %>%
  # CCPFoundingDay --> CCP_Founding_Day
  rename_with(~ str_replace(., "CCPFoundingDay$", "CCP_Founding_Day"),
            starts_with("t_")) %>%
  rename_with(~ str_replace(., "CCPFounding$", "CCP_Founding_Day"),
              starts_with("t_"))


# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

# Function 1: Extract odds ratios with standard (naive) SEs -------------------
extract_or_derived <- function(model, target_outcome, reference_outcome, conf_level = 0.95) {
  # Get coefficients and variance-covariance matrix
  coefs <- coef(model)
  vcov_mat <- vcov(model)
  
  # Get base category and all outcome names
  base_outcome <- model$lev[1]  # First level is the baseline
  outcome_names <- rownames(coefs)  # Non-baseline outcomes
  all_outcomes <- c(base_outcome, outcome_names)  # All possible outcomes
  
  # Get variable names (excluding intercept and fixed effects)
  var_names <- colnames(coefs)
  var_names <- var_names[!grepl("^\\(Intercept\\)|^city_name|^year", var_names)]
  
  # Check if target and reference outcomes exist in model
  if (!target_outcome %in% all_outcomes) {
    stop(paste("Target outcome", target_outcome, "not found in model. Available outcomes:", 
               paste(all_outcomes, collapse = ", ")))
  }
  if (!reference_outcome %in% all_outcomes) {
    stop(paste("Reference outcome", reference_outcome, "not found in model. Available outcomes:",
               paste(all_outcomes, collapse = ", ")))
  }
  
  # Calculate z critical value
  z <- qnorm((1 + conf_level) / 2)
  
  # Initialize results
  or_data <- tibble()
  
  # Loop through variables
  for (var_name in var_names) {
    var_idx <- which(colnames(coefs) == var_name)
    
    # Determine if target and reference are base category
    target_is_base <- (target_outcome == base_outcome)
    ref_is_base <- (reference_outcome == base_outcome)
    
    if (target_is_base && ref_is_base) {
      # Both are base - comparison is meaningless (0 log-odds)
      stop("Cannot compare base category to itself")
      
    } else if (target_is_base && !ref_is_base) {
      # Target is base, reference is not
      # log-odds(base vs ref) = -log-odds(ref vs base) = -beta_ref
      ref_idx <- which(rownames(coefs) == reference_outcome)
      beta_ref <- coefs[ref_idx, var_idx]
      
      beta_derived <- -beta_ref
      
      # Variance is the same, just SE of reference coefficient
      n_vars <- ncol(coefs)
      vcov_idx_ref <- (ref_idx - 1) * n_vars + var_idx
      se_derived <- sqrt(vcov_mat[vcov_idx_ref, vcov_idx_ref])
      
    } else if (!target_is_base && ref_is_base) {
      # Target is not base, reference is base
      # log-odds(target vs base) = beta_target (already relative to base)
      target_idx <- which(rownames(coefs) == target_outcome)
      beta_target <- coefs[target_idx, var_idx]
      
      beta_derived <- beta_target
      
      # Variance is just the variance of target coefficient
      n_vars <- ncol(coefs)
      vcov_idx_target <- (target_idx - 1) * n_vars + var_idx
      se_derived <- sqrt(vcov_mat[vcov_idx_target, vcov_idx_target])
      
    } else {
      # Both are non-base categories
      # log-odds(target vs ref) = beta_target - beta_ref
      target_idx <- which(rownames(coefs) == target_outcome)
      ref_idx <- which(rownames(coefs) == reference_outcome)
      
      beta_target <- coefs[target_idx, var_idx]
      beta_ref <- coefs[ref_idx, var_idx]
      
      beta_derived <- beta_target - beta_ref
      
      # Calculate variance of derived coefficient using vcov matrix
      # Var(β_target - β_ref) = Var(β_target) + Var(β_ref) - 2*Cov(β_target, β_ref)
      n_vars <- ncol(coefs)
      vcov_idx_target <- (target_idx - 1) * n_vars + var_idx
      vcov_idx_ref <- (ref_idx - 1) * n_vars + var_idx
      
      var_target <- vcov_mat[vcov_idx_target, vcov_idx_target]
      var_ref <- vcov_mat[vcov_idx_ref, vcov_idx_ref]
      cov_target_ref <- vcov_mat[vcov_idx_target, vcov_idx_ref]
      
      se_derived <- sqrt(var_target + var_ref - 2 * cov_target_ref)
    }
    
    # Calculate odds ratio and confidence intervals
    or <- exp(beta_derived)
    or_lower <- exp(beta_derived - z * se_derived)
    or_upper <- exp(beta_derived + z * se_derived)
    
    # Add to results
    or_data <- bind_rows(
      or_data,
      tibble(
        outcome = paste(target_outcome, "vs", reference_outcome),
        variable = var_name,
        or = or,
        conf_lower = or_lower,
        conf_upper = or_upper,
        conf_level = conf_level,
        se_type = "naive"
      )
    )
  }
  
  return(or_data)
}


# Function 2: Cluster-robust variance-covariance matrix -----------------------
vcovCL_multinom <- function(model, cluster) {
  
  cat("  Computing cluster-robust variance-covariance matrix...\n")
  
  # 1. Extract basic model components
  coefs <- coef(model)
  if (is.vector(coefs)) {
    coefs <- matrix(coefs, nrow = 1)
  }
  n_outcomes <- nrow(coefs)
  n_vars <- ncol(coefs)
  
  # 2. Get model matrix and model frame
  X <- model.matrix(model)
  mf <- model.frame(model)
  
  # 3. Get fitted probabilities
  fitted_probs <- fitted(model)
  if (is.vector(fitted_probs)) {
    fitted_probs <- matrix(fitted_probs, ncol = length(model$lev))
  }
  
  # 4. Extract actual response variable from model frame
  y_actual <- model.response(mf)
  if (is.factor(y_actual)) {
    y_actual <- as.character(y_actual)
  }
  
  y_levels <- model$lev
  
  # 5. Handle missing data: align cluster with model observations
  if (length(cluster) != nrow(X)) {
    if (!is.null(rownames(X))) {
      if (length(cluster) > nrow(X)) {
        cat("  Adjusting cluster for missing data...\n")
        cluster_indices <- as.numeric(rownames(X))
        cluster_vec <- cluster[cluster_indices]
      } else {
        stop("Cluster vector shorter than model data")
      }
    } else if (!is.null(model$na.action)) {
      cat("  Adjusting for missing data via na.action...\n")
      na_indices <- as.numeric(model$na.action)
      all_indices <- seq_along(cluster)
      kept_indices <- all_indices[-na_indices]
      cluster_vec <- cluster[kept_indices]
      cat("  Dropped", length(na_indices), "observations with missing values\n")
    } else {
      stop("Cannot align cluster vector with model data")
    }
  } else {
    cluster_vec <- cluster
  }
  
  # Verify alignment
  if (length(cluster_vec) != nrow(X)) {
    stop("Cluster vector length mismatch after adjustment")
  }
  
  # 6. Create indicator matrix for actual outcomes
  Y <- matrix(0, nrow = length(y_actual), ncol = length(y_levels))
  colnames(Y) <- y_levels
  
  for (i in seq_along(y_levels)) {
    Y[, i] <- as.numeric(y_actual == y_levels[i])
  }
  
  # 7. Calculate residuals
  residuals <- Y - fitted_probs
  
  # 8. Get standard variance-covariance matrix (bread) directly from model
  bread <- vcov(model)
  
  # Ensure bread is symmetric
  if (!isSymmetric(bread, tol = 1e-10)) {
    warning("Standard vcov not symmetric - symmetrizing")
    bread <- (bread + t(bread)) / 2
  }
  
  # 9. Calculate cluster-robust meat
  cluster_vec <- as.character(cluster_vec)
  cluster_ids <- unique(cluster_vec)
  n_clusters <- length(cluster_ids)
  
  cat("  Number of clusters:", n_clusters, "\n")
  
  if (n_clusters < 30) {
    warning("Fewer than 30 clusters - asymptotic approximation may be poor")
  }
  
  # Initialize meat matrix
  n_params <- n_outcomes * n_vars
  meat <- matrix(0, nrow = n_params, ncol = n_params)
  
  # Loop over clusters
  for (c in cluster_ids) {
    in_cluster <- which(cluster_vec == c)
    
    X_c <- X[in_cluster, , drop = FALSE]
    resid_c <- residuals[in_cluster, , drop = FALSE]
    
    # Build score vector for this cluster
    score_c <- numeric(n_params)
    
    for (j in 1:n_outcomes) {
      idx_start <- (j - 1) * n_vars + 1
      idx_end <- j * n_vars
      
      # Score = X' * residual for this outcome
      # Use j + 1 because first column is baseline (skipped)
      score_c[idx_start:idx_end] <- crossprod(X_c, resid_c[, j + 1])
    }
    
    # Add outer product to meat
    meat <- meat + tcrossprod(score_c)
  }
  
  # 10. Finite sample adjustment
  n_obs <- nrow(X)
  k <- n_params
  
  # Standard finite sample correction
  adjustment <- (n_clusters / (n_clusters - 1)) * ((n_obs - 1) / (n_obs - k))
  meat <- meat * adjustment
  
  # 11. Compute sandwich variance-covariance matrix
  vcov_cl <- bread %*% meat %*% bread
  
  # Ensure symmetry
  vcov_cl <- (vcov_cl + t(vcov_cl)) / 2
  
  # Check for numerical issues
  if (any(!is.finite(vcov_cl))) {
    warning("Non-finite values in cluster-robust vcov")
    vcov_cl[!is.finite(vcov_cl)] <- 0
  }
  
  cat("  Cluster-robust variance-covariance matrix computed.\n")
  return(vcov_cl)
}


# Function 3: Extract OR with cluster-robust SEs -------------------------------
extract_or_derived_clustered <- function(model, cluster, target_outcome, reference_outcome, 
                                         conf_level = 0.95) {
  
  # Calculate cluster-robust vcov
  vcov_cl <- vcovCL_multinom(model, cluster)
  
  # Get coefficients
  coefs <- coef(model)
  
  # Get base category and all outcome names
  base_outcome <- model$lev[1]
  outcome_names <- rownames(coefs)
  all_outcomes <- c(base_outcome, outcome_names)
  
  # Get variable names (excluding intercept and fixed effects)
  var_names <- colnames(coefs)
  var_names <- var_names[!grepl("^\\(Intercept\\)|^city_name|^year", var_names)]
  
  # Check if target and reference outcomes exist in model
  if (!target_outcome %in% all_outcomes) {
    stop(paste("Target outcome", target_outcome, "not found in model. Available outcomes:", 
               paste(all_outcomes, collapse = ", ")))
  }
  if (!reference_outcome %in% all_outcomes) {
    stop(paste("Reference outcome", reference_outcome, "not found in model. Available outcomes:",
               paste(all_outcomes, collapse = ", ")))
  }
  
  # Calculate z critical value
  z <- qnorm((1 + conf_level) / 2)
  
  # Initialize results
  or_data <- tibble()
  
  # Loop through variables
  for (var_name in var_names) {
    var_idx <- which(colnames(coefs) == var_name)
    
    # Determine if target and reference are base category
    target_is_base <- (target_outcome == base_outcome)
    ref_is_base <- (reference_outcome == base_outcome)
    
    if (target_is_base && ref_is_base) {
      stop("Cannot compare base category to itself")
      
    } else if (target_is_base && !ref_is_base) {
      # Target is base, reference is not
      ref_idx <- which(rownames(coefs) == reference_outcome)
      beta_ref <- coefs[ref_idx, var_idx]
      
      beta_derived <- -beta_ref
      
      # Variance from cluster-robust vcov
      n_vars <- ncol(coefs)
      vcov_idx_ref <- (ref_idx - 1) * n_vars + var_idx
      se_derived <- sqrt(vcov_cl[vcov_idx_ref, vcov_idx_ref])
      
    } else if (!target_is_base && ref_is_base) {
      # Target is not base, reference is base
      target_idx <- which(rownames(coefs) == target_outcome)
      beta_target <- coefs[target_idx, var_idx]
      
      beta_derived <- beta_target
      
      # Variance from cluster-robust vcov
      n_vars <- ncol(coefs)
      vcov_idx_target <- (target_idx - 1) * n_vars + var_idx
      se_derived <- sqrt(vcov_cl[vcov_idx_target, vcov_idx_target])
      
    } else {
      # Both are non-base categories
      target_idx <- which(rownames(coefs) == target_outcome)
      ref_idx <- which(rownames(coefs) == reference_outcome)
      
      beta_target <- coefs[target_idx, var_idx]
      beta_ref <- coefs[ref_idx, var_idx]
      
      beta_derived <- beta_target - beta_ref
      
      # Calculate variance using cluster-robust vcov
      n_vars <- ncol(coefs)
      vcov_idx_target <- (target_idx - 1) * n_vars + var_idx
      vcov_idx_ref <- (ref_idx - 1) * n_vars + var_idx
      
      var_target <- vcov_cl[vcov_idx_target, vcov_idx_target]
      var_ref <- vcov_cl[vcov_idx_ref, vcov_idx_ref]
      cov_target_ref <- vcov_cl[vcov_idx_target, vcov_idx_ref]
      
      se_derived <- sqrt(var_target + var_ref - 2 * cov_target_ref)
    }
    
    # Calculate odds ratio and confidence intervals
    or <- exp(beta_derived)
    or_lower <- exp(beta_derived - z * se_derived)
    or_upper <- exp(beta_derived + z * se_derived)
    
    # Add to results
    or_data <- bind_rows(
      or_data,
      tibble(
        outcome = paste(target_outcome, "vs", reference_outcome),
        variable = var_name,
        or = or,
        conf_lower = or_lower,
        conf_upper = or_upper,
        conf_level = conf_level,
        se_type = "clustered"
      )
    )
  }
  
  return(or_data)
}



# ==============================================================================
# DATA PREPARATION
# ==============================================================================

cat("\n===== Data Preparation =====\n\n")

# Recode city names ------------------------------------------------------------
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
    ),
    city_name = factor(city_name, levels = sort(unique(city_name)))
  )

# Prepare variables ------------------------------------------------------------
data <- data %>%
  mutate(
    # Convert categorical variables to factors
    sc_signal_4cat = factor(sc_signal_4cat, 
                            levels = c("Coercive", "Coercive + Administrative", "Administrative", 
                                       "No Signaling")),
    sc_signal_3cat = factor(sc_signal_3cat,
                            levels = c("Coercive", "Administrative", "No Signaling")),
    year = factor(year),
    create_date = as.Date(create_date),
    
    # Create cluster variables
    cluster_city = city_name,  # For Part II
    cluster_city_year = paste(city_name, year, sep = "_")  # For Part III
  )

cat("Data prepared successfully.\n")
cat("Cities:", nlevels(data$city_name), "\n")
cat("Years:", nlevels(data$year), "\n")
cat("City clusters:", length(unique(data$cluster_city)), "\n")
cat("City-year clusters:", length(unique(data$cluster_city_year)), "\n\n")


# ==============================================================================
# PART I: PREVALENCE ANALYSIS (UNCHANGED)
# ==============================================================================

cat("===== PART I: Prevalence Analysis =====\n\n")

# 1. Prevalence by city --------------------------------------------------------

# 4-category
prev_city_4cat <- data %>%
  count(city_name, sc_signal_4cat) %>%
  group_by(city_name) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

plot_city_4cat <- prev_city_4cat %>%
  ggplot(aes(x = city_name, y = proportion, fill = sc_signal_4cat)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Set2", name = "Signal Type") +
  labs(
    title = "State Capacity Signaling by City",
    subtitle = "4-Category Classification",
    x = "City",
    y = "Proportion"
  ) +
  # Guide fill in two rows
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "part1_prevalence_by_city_4cat.pdf"), 
       plot_city_4cat, width = 10, height = 5, device = "pdf")
cat("Saved: part1_prevalence_by_city_4cat.pdf\n")

# 3-category
prev_city_3cat <- data %>%
  count(city_name, sc_signal_3cat) %>%
  group_by(city_name) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup() %>%
  mutate(
    sc_signal_3cat = factor(sc_signal_3cat,
                            levels = c("Coercive", "No Signaling", "Administrative"))
  )

plot_city_3cat <- prev_city_3cat %>%
  ggplot(aes(x = city_name, y = proportion, fill = sc_signal_3cat)) +
  geom_col(position = "fill", alpha = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Set1", name = "Signal Type") +
  labs(
    title = "State Capacity Signaling by City",
    x = "City",
    y = "Proportion"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "part1_prevalence_by_city_3cat.pdf"), 
       plot_city_3cat, width = 10, height = 5, device = "pdf")
cat("Saved: part1_prevalence_by_city_3cat.pdf\n")


# 2. Prevalence over time ------------------------------------------------------

# 4-category
prev_time_4cat <- data %>%
  mutate(year = format(create_date, "%Y-%m-01")) %>%
  mutate(year = as.Date(year)) %>%
  count(year, sc_signal_4cat) %>%
  # count(year, sc_signal_4cat) %>%
  group_by(year) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

plot_time_4cat <- prev_time_4cat %>%
  ggplot(aes(x = year, y = proportion, color = sc_signal_4cat, 
             group = sc_signal_4cat)) +
  # geom_line(linewidth = 1.2) +
  # Add a smooth line for each category, no confidence interval
  geom_smooth(linetype = "dashed", alpha = 0.5) +
  geom_point(aes(shape = sc_signal_4cat), size = 2, alpha = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "State Capacity Signaling Over Time",
    subtitle = "4-Category Classification",
    x = "Year",
    y = "Proportion",
    color = "Signal Type", shape = "Signal Type"
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  guides(color = guide_legend(nrow=2, byrow=TRUE),
         shape = guide_legend(nrow=2, byrow=TRUE)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "part1_prevalence_by_time_4cat.pdf"), 
       plot_time_4cat, width = 8, height = 5, device = "pdf")
cat("Saved: part1_prevalence_by_time_4cat.pdf\n")

# 3-category
prev_time_3cat <- data %>%
  mutate(year = format(create_date, "%Y-%m-01")) %>%
  mutate(year = as.Date(year)) %>%
  count(year, sc_signal_3cat) %>%
  group_by(year) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup() %>%
  mutate(
    sc_signal_3cat = factor(sc_signal_3cat,
                            levels = c("Coercive", "No Signaling", "Administrative"))
  )

plot_time_3cat <- prev_time_3cat %>%
  ggplot(aes(x = year, y = proportion, color = sc_signal_3cat, 
             group = sc_signal_3cat)) +
  # geom_line(linewidth = 1.2) +
  # Add a smooth line for each category, no confidence interval
  geom_smooth(linetype = "dashed", alpha = 0.3) +
  geom_point(aes(shape = sc_signal_3cat), size = 2, alpha = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "State Capacity Signaling Over Time",
    subtitle = "3-Category Classification",
    x = "Year",
    y = "Proportion",
    color = "Signal Type", shape = "Signal Type"
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  guides(color = guide_legend(nrow=1, byrow=TRUE),
         shape = guide_legend(nrow=1, byrow=TRUE)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "part1_prevalence_by_time_3cat.pdf"), 
       plot_time_3cat, width = 8, height = 5, device = "pdf")
cat("Saved: part1_prevalence_by_time_3cat.pdf\n\n")


# ==============================================================================
# PART II: EVENT EFFECTS ANALYSIS (WITH CLUSTER-ROBUST SEs)
# ==============================================================================

cat("===== PART II: Event Effects Analysis =====\n\n")
cat("Clustering standard errors by CITY\n\n")

# Define event timing categories (reordered: 14d prior → during)
event_timings <- c("prior14d", "prior7d", "prior3d", "during")
timing_labels <- c("14 Days Prior", "7 Days Prior", "3 Days Prior", "During Event")

# Prepare Part II plot data only.
# This filters plots to the selected events, but leaves models/tables unchanged.

# Events to show in Part II plots only
# Desired order:
# Chinese New Year, National Day, CCP Founding Day, Party Congress, Two Sessions
part2_plot_events <- c(
  "Chinese_New_Year",
  "National_Day",
  "CCP_Founding_Day",
  "Party_Congress",
  "Two_Sessions"
)

part2_plot_event_labels <- c(
  "Chinese New Year",
  "National Day",
  "CCP Founding Day",
  "Party Congress",
  "Two Sessions"
)

part2_plot_event_lookup <- tibble(
  variable = part2_plot_events,
  event_clean = part2_plot_event_labels
)

prep_part2_plot_data <- function(x) {
  x %>%
    mutate(
      variable = str_remove(
        variable,
        "^t_(prior14d|prior7d|prior3d|during)_"
      ),
      
      # Important fix: harmonize CCP Founding naming
      variable = case_when(
        variable %in% c("CCPFounding", "CCPFoundingDay") ~ "CCP_Founding_Day",
        TRUE ~ variable
      ),
      
      timing = case_when(
        as.character(timing) == "prior14d" ~ "14 Days Prior",
        as.character(timing) == "prior7d"  ~ "7 Days Prior",
        as.character(timing) == "prior3d"  ~ "3 Days Prior",
        as.character(timing) == "during"   ~ "During Event",
        TRUE ~ as.character(timing)
      )
    ) %>%
    filter(variable %in% part2_plot_events) %>%
    left_join(part2_plot_event_lookup, by = "variable") %>%
    mutate(
      event_clean = factor(event_clean, levels = part2_plot_event_labels),
      timing = factor(timing, levels = timing_labels),
      ci_level = factor(ci_level, levels = c("90%", "95%")),
      se_type = factor(se_type, levels = c("naive", "clustered"))
    )
}



# Initialize storage for results
part2_models_3cat <- list()
part2_models_4cat <- list()

part2_results_3cat_naive_95 <- list()
part2_results_3cat_clustered_95 <- list()
part2_results_3cat_naive_90 <- list()
part2_results_3cat_clustered_90 <- list()

part2_results_4cat_admin_naive_95 <- list()
part2_results_4cat_admin_clustered_95 <- list()
part2_results_4cat_admin_naive_90 <- list()
part2_results_4cat_admin_clustered_90 <- list()

part2_results_4cat_both_naive_95 <- list()
part2_results_4cat_both_clustered_95 <- list()
part2_results_4cat_both_naive_90 <- list()
part2_results_4cat_both_clustered_90 <- list()

# Loop through each timing category --------------------------------------------
for (timing_name in event_timings) {
  
  cat("\n--- Analyzing:", timing_name, "---\n")
  
  # Select relevant event variables
  pattern <- paste0("^t_", timing_name, "_")
  event_vars <- data %>%
    select(matches(pattern)) %>%
    names()
  
  if (length(event_vars) == 0) {
    cat("No variables found for pattern:", pattern, "\n")
    next
  }
  
  cat("Variables:", length(event_vars), "\n")
  
  # Create renamed data for this timing
  # Remove timing prefix so texreg can group same events across models
  data_renamed <- data
  for (var in event_vars) {
    new_name <- str_remove(var, paste0("^t_", timing_name, "_"))
    data_renamed[[new_name]] <- data_renamed[[var]]
  }
  
  # Get renamed variable names
  event_vars_renamed <- str_remove(event_vars, paste0("^t_", timing_name, "_"))
  
  # Model formulas using renamed variables
  formula_3cat <- as.formula(paste(
    "sc_signal_3cat ~", 
    paste(event_vars_renamed, collapse = " + "),
    "+ city_name + year"
  ))
  
  formula_4cat <- as.formula(paste(
    "sc_signal_4cat ~", 
    paste(event_vars_renamed, collapse = " + "),
    "+ city_name + year"
  ))
  
  # Fit 3-category model
  cat("Fitting 3-category model...\n")
  model_3cat <- multinom(formula_3cat, data = data_renamed, 
                         trace = FALSE, model = TRUE)
  part2_models_3cat[[timing_name]] <- model_3cat
  
  # Fit 4-category model
  cat("Fitting 4-category model...\n")
  model_4cat <- multinom(formula_4cat, data = data_renamed, 
                         trace = FALSE, model = TRUE)
  part2_models_4cat[[timing_name]] <- model_4cat
  
  # Extract OR for 3-category model with 95% CI
  cat("Extracting 3-category odds ratios with naive SEs (95% CI)...\n")
  or_3cat_naive_95 <- extract_or_derived(
    model = model_3cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting 3-category odds ratios with cluster-robust SEs (95% CI)...\n")
  or_3cat_clustered_95 <- extract_or_derived_clustered(
    model = model_3cat,
    cluster = data_renamed$cluster_city,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  # Extract OR for 3-category model with 90% CI
  cat("Extracting 3-category odds ratios with naive SEs (90% CI)...\n")
  or_3cat_naive_90 <- extract_or_derived(
    model = model_3cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  cat("Extracting 3-category odds ratios with cluster-robust SEs (90% CI)...\n")
  or_3cat_clustered_90 <- extract_or_derived_clustered(
    model = model_3cat,
    cluster = data_renamed$cluster_city,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Extract OR for 4-category model: Administrative vs Coercive
  cat("Extracting 4-category odds ratios (Admin vs Coercive) with naive SEs (95% CI)...\n")
  or_4cat_admin_naive_95 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting 4-category odds ratios (Admin vs Coercive) with cluster-robust SEs (95% CI)...\n")
  or_4cat_admin_clustered_95 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_renamed$cluster_city,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting 4-category odds ratios (Admin vs Coercive) with naive SEs (90% CI)...\n")
  or_4cat_admin_naive_90 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  cat("Extracting 4-category odds ratios (Admin vs Coercive) with cluster-robust SEs (90% CI)...\n")
  or_4cat_admin_clustered_90 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_renamed$cluster_city,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Extract OR for 4-category model: Coercive+Administrative vs Coercive
  cat("Extracting 4-category odds ratios (Both vs Coercive) with naive SEs (95% CI)...\n")
  or_4cat_both_naive_95 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting 4-category odds ratios (Both vs Coercive) with cluster-robust SEs (95% CI)...\n")
  or_4cat_both_clustered_95 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_renamed$cluster_city,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting 4-category odds ratios (Both vs Coercive) with naive SEs (90% CI)...\n")
  or_4cat_both_naive_90 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  cat("Extracting 4-category odds ratios (Both vs Coercive) with cluster-robust SEs (90% CI)...\n")
  or_4cat_both_clustered_90 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_renamed$cluster_city,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Store results
  part2_results_3cat_naive_95[[timing_name]] <- or_3cat_naive_95 %>% mutate(timing = timing_name)
  part2_results_3cat_clustered_95[[timing_name]] <- or_3cat_clustered_95 %>% mutate(timing = timing_name)
  part2_results_3cat_naive_90[[timing_name]] <- or_3cat_naive_90 %>% mutate(timing = timing_name)
  part2_results_3cat_clustered_90[[timing_name]] <- or_3cat_clustered_90 %>% mutate(timing = timing_name)
  
  part2_results_4cat_admin_naive_95[[timing_name]] <- or_4cat_admin_naive_95 %>% mutate(timing = timing_name)
  part2_results_4cat_admin_clustered_95[[timing_name]] <- or_4cat_admin_clustered_95 %>% mutate(timing = timing_name)
  part2_results_4cat_admin_naive_90[[timing_name]] <- or_4cat_admin_naive_90 %>% mutate(timing = timing_name)
  part2_results_4cat_admin_clustered_90[[timing_name]] <- or_4cat_admin_clustered_90 %>% mutate(timing = timing_name)
  
  part2_results_4cat_both_naive_95[[timing_name]] <- or_4cat_both_naive_95 %>% mutate(timing = timing_name)
  part2_results_4cat_both_clustered_95[[timing_name]] <- or_4cat_both_clustered_95 %>% mutate(timing = timing_name)
  part2_results_4cat_both_naive_90[[timing_name]] <- or_4cat_both_naive_90 %>% mutate(timing = timing_name)
  part2_results_4cat_both_clustered_90[[timing_name]] <- or_4cat_both_clustered_90 %>% mutate(timing = timing_name)
}

# Combine all results ----------------------------------------------------------
part2_3cat_naive_95 <- bind_rows(part2_results_3cat_naive_95)
part2_3cat_clustered_95 <- bind_rows(part2_results_3cat_clustered_95)
part2_3cat_naive_90 <- bind_rows(part2_results_3cat_naive_90)
part2_3cat_clustered_90 <- bind_rows(part2_results_3cat_clustered_90)

part2_4cat_admin_naive_95 <- bind_rows(part2_results_4cat_admin_naive_95)
part2_4cat_admin_clustered_95 <- bind_rows(part2_results_4cat_admin_clustered_95)
part2_4cat_admin_naive_90 <- bind_rows(part2_results_4cat_admin_naive_90)
part2_4cat_admin_clustered_90 <- bind_rows(part2_results_4cat_admin_clustered_90)

part2_4cat_both_naive_95 <- bind_rows(part2_results_4cat_both_naive_95)
part2_4cat_both_clustered_95 <- bind_rows(part2_results_4cat_both_clustered_95)
part2_4cat_both_naive_90 <- bind_rows(part2_results_4cat_both_naive_90)
part2_4cat_both_clustered_90 <- bind_rows(part2_results_4cat_both_clustered_90)

# Combine all for plotting (both CI levels)
part2_3cat_combined_all <- bind_rows(
  part2_3cat_naive_95 %>% mutate(ci_level = "95%"),
  part2_3cat_clustered_95 %>% mutate(ci_level = "95%"),
  part2_3cat_naive_90 %>% mutate(ci_level = "90%"),
  part2_3cat_clustered_90 %>% mutate(ci_level = "90%")
)

part2_4cat_admin_combined_all <- bind_rows(
  part2_4cat_admin_naive_95 %>% mutate(ci_level = "95%"),
  part2_4cat_admin_clustered_95 %>% mutate(ci_level = "95%"),
  part2_4cat_admin_naive_90 %>% mutate(ci_level = "90%"),
  part2_4cat_admin_clustered_90 %>% mutate(ci_level = "90%")
)

part2_4cat_both_combined_all <- bind_rows(
  part2_4cat_both_naive_95 %>% mutate(ci_level = "95%"),
  part2_4cat_both_clustered_95 %>% mutate(ci_level = "95%"),
  part2_4cat_both_naive_90 %>% mutate(ci_level = "90%"),
  part2_4cat_both_clustered_90 %>% mutate(ci_level = "90%")
)

part2_3cat_plot_all <- prep_part2_plot_data(part2_3cat_combined_all)

part2_4cat_admin_plot_all <- prep_part2_plot_data(part2_4cat_admin_combined_all)

part2_4cat_both_plot_all <- prep_part2_plot_data(part2_4cat_both_combined_all)


# Create LaTeX tables ----------------------------------------------------------
cat("\n--- Creating LaTeX tables ---\n")

# Get coefficient names from first model
all_coef_names_3cat <- names(coef(part2_models_3cat[[1]]))
custom_coef_names_3cat <- all_coef_names_3cat
names(custom_coef_names_3cat) <- all_coef_names_3cat

# Clean up event names (they no longer have timing prefixes)
for (i in seq_along(all_coef_names_3cat)) {
  if (!grepl("^\\(Intercept\\)|^city_name|^year", all_coef_names_3cat[i])) {
    custom_coef_names_3cat[i] <- all_coef_names_3cat[i] %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  }
}

# Table 1: 3-category models with naive SEs
texreg(
  part2_models_3cat,
  file = file.path(output_dir, "part2_3cat_models_naive_se.tex"),
  custom.model.names = timing_labels,
  custom.coef.names = custom_coef_names_3cat,
  omit.coef = "city_name|year",
  caption = "Event Effects on State Capacity Signaling (Naive Standard Errors)",
  label = "tab:part2_3cat_naive",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = "Standard errors in parentheses. %stars. Administrative vs Coercive signaling.",
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "City FE" = rep("Yes", length(part2_models_3cat)),
    "Year FE" = rep("Yes", length(part2_models_3cat))
  )
)
cat("Saved: part2_3cat_models_naive_se.tex\n")

# Table 2: 3-category models with cluster-robust SEs (WITH CORRECTED STARS)
clustered_ses_3cat <- lapply(names(part2_models_3cat), function(timing_name) {
  model <- part2_models_3cat[[timing_name]]
  pattern <- paste0("^t_", timing_name, "_")
  event_vars <- data %>% select(matches(pattern)) %>% names()
  data_renamed <- data
  for (var in event_vars) {
    new_name <- str_remove(var, paste0("^t_", timing_name, "_"))
    data_renamed[[new_name]] <- data_renamed[[var]]
  }
  vcov_cl <- vcovCL_multinom(model, data_renamed$cluster_city)
  sqrt(diag(vcov_cl))
})

# Calculate corrected p-values for clustered SEs
clustered_pvalues_3cat <- lapply(names(part2_models_3cat), function(timing_name) {
  model <- part2_models_3cat[[timing_name]]
  pattern <- paste0("^t_", timing_name, "_")
  event_vars <- data %>% select(matches(pattern)) %>% names()
  data_renamed <- data
  for (var in event_vars) {
    new_name <- str_remove(var, paste0("^t_", timing_name, "_"))
    data_renamed[[new_name]] <- data_renamed[[var]]
  }
  vcov_cl <- vcovCL_multinom(model, data_renamed$cluster_city)
  clustered_se <- sqrt(diag(vcov_cl))
  
  # Get coefficients as vector
  coef_vec <- as.vector(t(coef(model)))
  
  # Calculate z-statistics and p-values
  z_stat <- coef_vec / clustered_se
  p_val <- 2 * pnorm(-abs(z_stat))
  
  return(p_val)
})

texreg(
  part2_models_3cat,
  override.se = clustered_ses_3cat,
  override.pvalues = clustered_pvalues_3cat,
  file = file.path(output_dir, "part2_3cat_models_clustered_se.tex"),
  custom.model.names = timing_labels,
  custom.coef.names = custom_coef_names_3cat,
  omit.coef = "city_name|year",
  caption = "Event Effects on State Capacity Signaling (Cluster-Robust Standard Errors)",
  label = "tab:part2_3cat_clustered",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = "Cluster-robust standard errors (by city) in parentheses. %stars. Administrative vs Coercive signaling.",
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "City FE" = rep("Yes", length(part2_models_3cat)),
    "Year FE" = rep("Yes", length(part2_models_3cat))
  )
)
cat("Saved: part2_3cat_models_clustered_se.tex\n")

# Tables for 4-category models
all_coef_names_4cat <- names(coef(part2_models_4cat[[1]]))
custom_coef_names_4cat <- all_coef_names_4cat
names(custom_coef_names_4cat) <- all_coef_names_4cat

for (i in seq_along(all_coef_names_4cat)) {
  if (!grepl("^\\(Intercept\\)|^city_name|^year", all_coef_names_4cat[i])) {
    custom_coef_names_4cat[i] <- all_coef_names_4cat[i] %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  }
}

# Table 3: 4-category models with naive SEs
texreg(
  part2_models_4cat,
  file = file.path(output_dir, "part2_4cat_models_naive_se.tex"),
  custom.model.names = timing_labels,
  custom.coef.names = custom_coef_names_4cat,
  omit.coef = "city_name|year",
  caption = "Event Effects on State Capacity Signaling - 4-Category Outcome (Naive Standard Errors)",
  label = "tab:part2_4cat_naive",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = "Standard errors in parentheses. %stars. 4-category outcome variable with separate Coercive+Administrative category.",
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "City FE" = rep("Yes", length(part2_models_4cat)),
    "Year FE" = rep("Yes", length(part2_models_4cat))
  )
)
cat("Saved: part2_4cat_models_naive_se.tex\n")

# Table 4: 4-category models with cluster-robust SEs (WITH CORRECTED STARS)
clustered_ses_4cat <- lapply(names(part2_models_4cat), function(timing_name) {
  model <- part2_models_4cat[[timing_name]]
  pattern <- paste0("^t_", timing_name, "_")
  event_vars <- data %>% select(matches(pattern)) %>% names()
  data_renamed <- data
  for (var in event_vars) {
    new_name <- str_remove(var, paste0("^t_", timing_name, "_"))
    data_renamed[[new_name]] <- data_renamed[[var]]
  }
  vcov_cl <- vcovCL_multinom(model, data_renamed$cluster_city)
  sqrt(diag(vcov_cl))
})

clustered_pvalues_4cat <- lapply(names(part2_models_4cat), function(timing_name) {
  model <- part2_models_4cat[[timing_name]]
  pattern <- paste0("^t_", timing_name, "_")
  event_vars <- data %>% select(matches(pattern)) %>% names()
  data_renamed <- data
  for (var in event_vars) {
    new_name <- str_remove(var, paste0("^t_", timing_name, "_"))
    data_renamed[[new_name]] <- data_renamed[[var]]
  }
  vcov_cl <- vcovCL_multinom(model, data_renamed$cluster_city)
  clustered_se <- sqrt(diag(vcov_cl))
  
  coef_vec <- as.vector(t(coef(model)))
  z_stat <- coef_vec / clustered_se
  p_val <- 2 * pnorm(-abs(z_stat))
  
  return(p_val)
})

texreg(
  part2_models_4cat,
  override.se = clustered_ses_4cat,
  override.pvalues = clustered_pvalues_4cat,
  file = file.path(output_dir, "part2_4cat_models_clustered_se.tex"),
  custom.model.names = timing_labels,
  custom.coef.names = custom_coef_names_4cat,
  omit.coef = "city_name|year",
  caption = "Event Effects on State Capacity Signaling - 4-Category Outcome (Cluster-Robust Standard Errors)",
  label = "tab:part2_4cat_clustered",
  single.row = FALSE,
  stars = c(0.01, 0.05, 0.1),
  custom.note = "Cluster-robust standard errors (by city) in parentheses. %stars. 4-category outcome variable with separate Coercive+Administrative category.",
  booktabs = TRUE,
  dcolumn = TRUE,
  use.packages = FALSE,
  custom.gof.rows = list(
    "City FE" = rep("Yes", length(part2_models_4cat)),
    "Year FE" = rep("Yes", length(part2_models_4cat))
  )
)
cat("Saved: part2_4cat_models_clustered_se.tex\n")


# Create unified plots ---------------------------------------------------------
cat("\n--- Creating unified plots ---\n")

# Helper function to create unified plots
create_unified_plot <- function(plot_data, se_type_label, title_suffix, 
                                output_filename, subtitle_text) {
  
  plot_data_prep <- plot_data %>%
    filter(se_type == se_type_label) %>%
    mutate(
      event_clean = factor(event_clean, levels = part2_plot_event_labels),
      timing = factor(timing, levels = timing_labels),
      ci_level = factor(ci_level, levels = c("90%", "95%"))
    )
  
  # Optional diagnostic: rows that cannot be shown on log scale
  bad_rows <- plot_data_prep %>%
    filter(
      !is.finite(or) | !is.finite(conf_lower) | !is.finite(conf_upper) |
        or <= 0 | conf_lower <= 0 | conf_upper <= 0
    )
  
  if (nrow(bad_rows) > 0) {
    warning("Some rows have non-positive or non-finite OR/CI values and may not plot on log scale.")
    print(bad_rows)
  }
  
  plot_data_95 <- plot_data_prep %>% filter(ci_level == "95%")
  plot_data_90 <- plot_data_prep %>% filter(ci_level == "90%")
  
  dodge <- position_dodge(width = 0.65)
  
  p <- ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    
    # 90% CI - thicker linerange
    geom_linerange(
      data = plot_data_90,
      aes(
        x = event_clean,
        ymin = conf_lower,
        ymax = conf_upper,
        color = timing,
        group = timing
      ),
      linewidth = 1,
      position = dodge
    ) +
    
    # 95% CI - thinner linerange
    geom_linerange(
      data = plot_data_95,
      aes(
        x = event_clean,
        ymin = conf_lower,
        ymax = conf_upper,
        color = timing,
        group = timing
      ),
      linewidth = 0.5,
      position = dodge
    ) +
    
    # Points
    geom_point(
      data = plot_data_95,
      aes(
        x = event_clean,
        y = or,
        color = timing,
        shape = timing,
        group = timing
      ),
      size = 3,
      position = dodge
    ) +
    
    scale_color_brewer(palette = "Set1", name = "Event Timing") +
    scale_shape_manual(values = c(16, 17, 15, 18), name = "Event Timing") +
    scale_x_discrete(drop = FALSE) +
    scale_y_log10() +
    labs(
      title = paste("Event Effects on State Capacity Signaling", title_suffix),
      subtitle = subtitle_text,
      x = "Event Type",
      y = "Odds Ratio",
      caption = ifelse(
        se_type_label == "naive",
        "Thick lines: 90% CI; Thin lines: 95% CI",
        "Thick lines: 90% CI; Thin lines: 95% CI; Standard errors clustered by city; Vertical axis on log scale"
      )
    ) +
    theme(
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom"
    )
  
  ggsave(
    file.path(output_dir, output_filename),
    p,
    width = 9,
    height = 6,
    device = "pdf"
  )
  
  cat("Saved:", output_filename, "\n")
  
  return(p)
}

# Plot 1: 3-category naive
create_unified_plot(
  plot_data = part2_3cat_plot_all,
  se_type_label = "naive",
  title_suffix = "",
  output_filename = "part2_3cat_unified_naive.pdf",
  subtitle_text = "Administrative vs Coercive Signaling (Naive SEs)"
)

# Plot 2: 3-category clustered
create_unified_plot(
  plot_data = part2_3cat_plot_all,
  se_type_label = "clustered",
  title_suffix = "",
  output_filename = "part2_3cat_unified_clustered.pdf",
  subtitle_text = "Administrative vs Coercive Signaling (Cluster-Robust SEs)"
)

# Plot 3: 4-category Admin vs Coercive naive
create_unified_plot(
  plot_data = part2_4cat_admin_plot_all,
  se_type_label = "naive",
  title_suffix = " - 4-Category",
  output_filename = "part2_4cat_admin_unified_naive.pdf",
  subtitle_text = "Administrative vs Coercive Signaling (Naive SEs)"
)

# Plot 4: 4-category Admin vs Coercive clustered
create_unified_plot(
  plot_data = part2_4cat_admin_plot_all,
  se_type_label = "clustered",
  title_suffix = " - 4-Category",
  output_filename = "part2_4cat_admin_unified_clustered.pdf",
  subtitle_text = "Administrative vs Coercive Signaling (Cluster-Robust SEs)"
)

# Plot 5: 4-category Both vs Coercive naive
create_unified_plot(
  plot_data = part2_4cat_both_plot_all,
  se_type_label = "naive",
  title_suffix = " - 4-Category",
  output_filename = "part2_4cat_both_unified_naive.pdf",
  subtitle_text = "Coercive+Administrative vs Coercive Signaling (Naive SEs)"
)

# Plot 6: 4-category Both vs Coercive clustered
create_unified_plot(
  plot_data = part2_4cat_both_plot_all,
  se_type_label = "clustered",
  title_suffix = " - 4-Category",
  output_filename = "part2_4cat_both_unified_clustered.pdf",
  subtitle_text = "Coercive+Administrative vs Coercive Signaling (Cluster-Robust SEs)"
)


# Save CSV files ---------------------------------------------------------------
write_csv(part2_3cat_naive_95, file.path(output_dir, "part2_3cat_naive_95ci.csv"))
write_csv(part2_3cat_clustered_95, file.path(output_dir, "part2_3cat_clustered_95ci.csv"))
write_csv(part2_3cat_naive_90, file.path(output_dir, "part2_3cat_naive_90ci.csv"))
write_csv(part2_3cat_clustered_90, file.path(output_dir, "part2_3cat_clustered_90ci.csv"))

write_csv(part2_4cat_admin_naive_95, file.path(output_dir, "part2_4cat_admin_naive_95ci.csv"))
write_csv(part2_4cat_admin_clustered_95, file.path(output_dir, "part2_4cat_admin_clustered_95ci.csv"))
write_csv(part2_4cat_admin_naive_90, file.path(output_dir, "part2_4cat_admin_naive_90ci.csv"))
write_csv(part2_4cat_admin_clustered_90, file.path(output_dir, "part2_4cat_admin_clustered_90ci.csv"))

write_csv(part2_4cat_both_naive_95, file.path(output_dir, "part2_4cat_both_naive_95ci.csv"))
write_csv(part2_4cat_both_clustered_95, file.path(output_dir, "part2_4cat_both_clustered_95ci.csv"))
write_csv(part2_4cat_both_naive_90, file.path(output_dir, "part2_4cat_both_naive_90ci.csv"))
write_csv(part2_4cat_both_clustered_90, file.path(output_dir, "part2_4cat_both_clustered_90ci.csv"))

cat("\nPart II analysis complete.\n\n")


# ==============================================================================
# PART III: CITY-LEVEL COVARIATES (WITH CLUSTER-ROBUST SEs)
# ==============================================================================

cat("===== PART III: City-Level Covariates Analysis =====\n\n")
cat("Clustering standard errors by CITY-YEAR\n\n")

# Create proportion variables (done once) --------------------------------------
data <- data %>%
  mutate(
    city_lag1_public_exp_publicservice_prop = 
      city_lag1_public_exp_publicservice_100000000yuan / city_lag1_public_exp_total_100000000yuan,
    city_lag1_public_exp_pupblicsecurity_prop = 
      city_lag1_public_exp_pupblicsecurity_100000000yuan / city_lag1_public_exp_total_100000000yuan
  )

# Helper function to clean variable names --------------------------------------
clean_var_names <- function(vars) {
  vars %>%
    str_remove("_scaled$") %>%
    str_remove("city_lag1_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

# Main function to fit multinomial models --------------------------------------
# ==============================================================================
# PART III: CITY-LEVEL COVARIATES (WITH CLUSTER-ROBUST SEs) - MODULAR VERSION
# ==============================================================================

cat("===== PART III: City-Level Covariates Analysis =====\n\n")
cat("Three specifications: Naive SEs, Clustered by CITY-YEAR, Clustered by CITY\n\n")

# Create proportion variables (done once) --------------------------------------
data <- data %>%
  mutate(
    city_lag1_public_exp_publicservice_prop = 
      city_lag1_public_exp_publicservice_100000000yuan / city_lag1_public_exp_total_100000000yuan,
    city_lag1_public_exp_pupblicsecurity_prop = 
      city_lag1_public_exp_pupblicsecurity_100000000yuan / city_lag1_public_exp_total_100000000yuan
  )

# Helper function to clean variable names --------------------------------------
clean_var_names <- function(vars) {
  vars %>%
    str_remove("_scaled$") %>%
    str_remove("city_lag1_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

# Main function to fit multinomial models --------------------------------------
fit_multinomial_models <- function(predictors, filename_suffix, title_label, 
                                   predictor_labels = predictor_labels_lookup) {
  
  cat("\n===== Fitting models:", title_label, "=====\n\n")
  cat("Predictors included:", length(predictors), "\n")
  cat(paste(predictors, collapse = "\n"), "\n\n")
  
  # Standardize variables ------------------------------------------------------
  data_scaled <- data %>%
    mutate(across(all_of(predictors), 
                  ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),
                  .names = "{.col}_scaled"))
  
  # Create mapping from scaled names to clean labels
  predictors_scaled <- paste0(predictors, "_scaled")
  
  rename_mapping <- setNames(
    predictors_scaled,  # old names (values)
    map_chr(predictors, function(pred) {  # new names (names of vector)
      label <- predictor_labels$clean_label[predictor_labels$predictor_name == pred]
      if (length(label) > 0) {
        return(label)
      } else {
        # Fallback if not in lookup table
        return(pred)
      }
    })
  )
  
  # Rename columns in data_scaled to use clean labels
  data_scaled <- data_scaled %>%
    rename(!!!rename_mapping)
  
  # Get the new predictor names (clean labels)
  predictors_labeled <- names(rename_mapping)
  
  # Model formulas using backtick-wrapped labels -----------------------------
  predictors_formula <- paste0("`", predictors_labeled, "`", collapse = " + ")
  
  formula_3cat <- as.formula(paste(
    "sc_signal_3cat ~", 
    predictors_formula,
    "+ year"
  ))
  
  formula_4cat <- as.formula(paste(
    "sc_signal_4cat ~", 
    predictors_formula,
    "+ year"
  ))
  
  # Fit models -----------------------------------------------------------------
  cat("Fitting 3-category model...\n")
  model_3cat <- multinom(formula_3cat, data = data_scaled, trace = FALSE)
  
  cat("Fitting 4-category model...\n")
  model_4cat <- multinom(formula_4cat, data = data_scaled, trace = FALSE)
  
  # Extract odds ratios for 3-category model ----------------------------------
  cat("\n--- 3-Category Model ---\n")
  
  # Naive SEs
  cat("Extracting odds ratios with naive SEs (95% CI)...\n")
  or_3cat_naive_95 <- extract_or_derived(
    model = model_3cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios with naive SEs (90% CI)...\n")
  or_3cat_naive_90 <- extract_or_derived(
    model = model_3cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Cluster-robust SEs (by City-Year)
  cat("Extracting odds ratios with cluster-robust SEs by City-Year (95% CI)...\n")
  or_3cat_clusteredCY_95 <- extract_or_derived_clustered(
    model = model_3cat,
    cluster = data_scaled$cluster_city_year,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios with cluster-robust SEs by City-Year (90% CI)...\n")
  or_3cat_clusteredCY_90 <- extract_or_derived_clustered(
    model = model_3cat,
    cluster = data_scaled$cluster_city_year,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Cluster-robust SEs (by City)
  cat("Extracting odds ratios with cluster-robust SEs by City (95% CI)...\n")
  or_3cat_clusteredC_95 <- extract_or_derived_clustered(
    model = model_3cat,
    cluster = data_scaled$city_name,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios with cluster-robust SEs by City (90% CI)...\n")
  or_3cat_clusteredC_90 <- extract_or_derived_clustered(
    model = model_3cat,
    cluster = data_scaled$city_name,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Extract odds ratios for 4-category model ----------------------------------
  cat("\n--- 4-Category Model ---\n")
  
  # Administrative vs Coercive - Naive SEs
  cat("Extracting odds ratios (Admin vs Coercive) with naive SEs (95% CI)...\n")
  or_4cat_admin_naive_95 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios (Admin vs Coercive) with naive SEs (90% CI)...\n")
  or_4cat_admin_naive_90 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Administrative vs Coercive - Cluster-robust SEs (by City-Year)
  cat("Extracting odds ratios (Admin vs Coercive) with cluster-robust SEs by City-Year (95% CI)...\n")
  or_4cat_admin_clusteredCY_95 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$cluster_city_year,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios (Admin vs Coercive) with cluster-robust SEs by City-Year (90% CI)...\n")
  or_4cat_admin_clusteredCY_90 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$cluster_city_year,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Administrative vs Coercive - Cluster-robust SEs (by City)
  cat("Extracting odds ratios (Admin vs Coercive) with cluster-robust SEs by City (95% CI)...\n")
  or_4cat_admin_clusteredC_95 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$city_name,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios (Admin vs Coercive) with cluster-robust SEs by City (90% CI)...\n")
  or_4cat_admin_clusteredC_90 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$city_name,
    target_outcome = "Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Coercive+Administrative vs Coercive - Naive SEs
  cat("Extracting odds ratios (Both vs Coercive) with naive SEs (95% CI)...\n")
  or_4cat_both_naive_95 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios (Both vs Coercive) with naive SEs (90% CI)...\n")
  or_4cat_both_naive_90 <- extract_or_derived(
    model = model_4cat,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Coercive+Administrative vs Coercive - Cluster-robust SEs (by City-Year)
  cat("Extracting odds ratios (Both vs Coercive) with cluster-robust SEs by City-Year (95% CI)...\n")
  or_4cat_both_clusteredCY_95 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$cluster_city_year,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios (Both vs Coercive) with cluster-robust SEs by City-Year (90% CI)...\n")
  or_4cat_both_clusteredCY_90 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$cluster_city_year,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Coercive+Administrative vs Coercive - Cluster-robust SEs (by City)
  cat("Extracting odds ratios (Both vs Coercive) with cluster-robust SEs by City (95% CI)...\n")
  or_4cat_both_clusteredC_95 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$city_name,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.95
  )
  
  cat("Extracting odds ratios (Both vs Coercive) with cluster-robust SEs by City (90% CI)...\n")
  or_4cat_both_clusteredC_90 <- extract_or_derived_clustered(
    model = model_4cat,
    cluster = data_scaled$city_name,
    target_outcome = "Coercive + Administrative",
    reference_outcome = "Coercive",
    conf_level = 0.90
  )
  
  # Create LaTeX tables --------------------------------------------------------
  cat("\n--- Creating LaTeX tables ---\n")
  
  # Get naive SEs and p-values for 3-category model
  cat("Computing all SE types for 3-category model...\n")
  naive_vcov_3cat <- vcov(model_3cat)
  naive_ses_3cat <- sqrt(diag(naive_vcov_3cat))
  coef_vec_3cat <- as.vector(t(coef(model_3cat)))
  z_stat_naive_3cat <- coef_vec_3cat / naive_ses_3cat
  naive_pvalues_3cat <- 2 * pnorm(-abs(z_stat_naive_3cat))
  
  # Get cluster-robust SEs for 3-category model
  vcov_cl_3cat_CY <- vcovCL_multinom(model_3cat, data_scaled$cluster_city_year)
  clustered_ses_3cat_CY <- sqrt(diag(vcov_cl_3cat_CY))
  z_stat_3cat_CY <- coef_vec_3cat / clustered_ses_3cat_CY
  clustered_pvalues_3cat_CY <- 2 * pnorm(-abs(z_stat_3cat_CY))
  
  vcov_cl_3cat_C <- vcovCL_multinom(model_3cat, data_scaled$city_name)
  clustered_ses_3cat_C <- sqrt(diag(vcov_cl_3cat_C))
  z_stat_3cat_C <- coef_vec_3cat / clustered_ses_3cat_C
  clustered_pvalues_3cat_C <- 2 * pnorm(-abs(z_stat_3cat_C))
  
  # Diagnostic output
  cat("\n  Diagnostic - lengths of SE vectors:\n")
  cat("    Naive:", length(naive_ses_3cat), "\n")
  cat("    Clustered (CY):", length(clustered_ses_3cat_CY), "\n")
  cat("    Clustered (C):", length(clustered_ses_3cat_C), "\n")
  cat("    Coefficients:", length(coef_vec_3cat), "\n")
  
  # Create custom coefficient names using lookup table
  all_coef_names <- colnames(coef(model_3cat))
  custom_coef_names <- setNames(all_coef_names, all_coef_names)
  
  for (i in seq_along(all_coef_names)) {
    coef_name <- all_coef_names[i]
    # Check if it's a predictor (remove _scaled suffix)
    pred_name <- str_remove(coef_name, "_scaled$")
    
    # Look up in the predictor labels table
    if (pred_name %in% predictor_labels$predictor_name) {
      label <- predictor_labels$clean_label[predictor_labels$predictor_name == pred_name]
      custom_coef_names[i] <- label
    }
  }
  
  # Three-column table: naive, clustered by City-Year, clustered by City
  texreg(
    list(model_3cat, model_3cat, model_3cat),
    override.se = list(naive_ses_3cat, clustered_ses_3cat_CY, clustered_ses_3cat_C),
    override.pvalues = list(naive_pvalues_3cat, clustered_pvalues_3cat_CY, clustered_pvalues_3cat_C),
    file = file.path(output_dir, paste0("part3_3cat_model_table_", filename_suffix, ".tex")),
    custom.model.names = c("Naive SE", "Clustered SE (by City-Year)", "Clustered SE (by City)"),
    # custom.coef.names = custom_coef_names,
    omit.coef = "year",
    caption = paste0("City-Level Determinants of State Capacity Signaling (", title_label, ")"),
    label = paste0("tab:part3_3cat_", filename_suffix),
    single.row = FALSE,
    stars = c(0.01, 0.05, 0.1),
    custom.note = "Standard errors in parentheses. %stars. Administrative vs Coercive signaling.",
    booktabs = TRUE,
    dcolumn = TRUE,
    use.packages = FALSE,
    custom.gof.rows = list(
      "Year FE" = c("Yes", "Yes", "Yes")
    )
  )
  cat("Saved: part3_3cat_model_table_", filename_suffix, ".tex\n", sep = "")
  
  # Get naive SEs and p-values for 4-category model
  cat("\nComputing all SE types for 4-category model...\n")
  naive_vcov_4cat <- vcov(model_4cat)
  naive_ses_4cat <- sqrt(diag(naive_vcov_4cat))
  coef_vec_4cat <- as.vector(t(coef(model_4cat)))
  z_stat_naive_4cat <- coef_vec_4cat / naive_ses_4cat
  naive_pvalues_4cat <- 2 * pnorm(-abs(z_stat_naive_4cat))
  
  # Get cluster-robust SEs for 4-category model
  vcov_cl_4cat_CY <- vcovCL_multinom(model_4cat, data_scaled$cluster_city_year)
  clustered_ses_4cat_CY <- sqrt(diag(vcov_cl_4cat_CY))
  z_stat_4cat_CY <- coef_vec_4cat / clustered_ses_4cat_CY
  clustered_pvalues_4cat_CY <- 2 * pnorm(-abs(z_stat_4cat_CY))
  
  vcov_cl_4cat_C <- vcovCL_multinom(model_4cat, data_scaled$city_name)
  clustered_ses_4cat_C <- sqrt(diag(vcov_cl_4cat_C))
  z_stat_4cat_C <- coef_vec_4cat / clustered_ses_4cat_C
  clustered_pvalues_4cat_C <- 2 * pnorm(-abs(z_stat_4cat_C))
  
  # Create custom coefficient names for 4-category model
  all_coef_names_4cat <- colnames(coef(model_4cat))
  custom_coef_names_4cat <- setNames(all_coef_names_4cat, all_coef_names_4cat)
  
  for (i in seq_along(all_coef_names_4cat)) {
    coef_name <- all_coef_names_4cat[i]
    pred_name <- str_remove(coef_name, "_scaled$")
    
    if (pred_name %in% predictor_labels$predictor_name) {
      label <- predictor_labels$clean_label[predictor_labels$predictor_name == pred_name]
      custom_coef_names_4cat[i] <- label
    }
  }
  
  # Three-column table: naive, clustered by City-Year, clustered by City
  texreg(
    list(model_4cat, model_4cat, model_4cat),
    override.se = list(naive_ses_4cat, clustered_ses_4cat_CY, clustered_ses_4cat_C),
    override.pvalues = list(naive_pvalues_4cat, clustered_pvalues_4cat_CY, clustered_pvalues_4cat_C),
    file = file.path(output_dir, paste0("part3_4cat_model_table_", filename_suffix, ".tex")),
    custom.model.names = c("Naive SE", "Clustered SE (by City-Year)", "Clustered SE (by City)"),
    # custom.coef.names = custom_coef_names_4cat,
    omit.coef = "year",
    caption = paste0("City-Level Determinants of State Capacity Signaling - 4-Category Outcome (", title_label, ")"),
    label = paste0("tab:part3_4cat_", filename_suffix),
    single.row = FALSE,
    stars = c(0.01, 0.05, 0.1),
    custom.note = "Standard errors in parentheses. %stars. 4-category outcome variable with separate Coercive+Administrative category.",
    booktabs = TRUE,
    dcolumn = TRUE,
    use.packages = FALSE,
    custom.gof.rows = list(
      "Year FE" = c("Yes", "Yes", "Yes")
    )
  )
  cat("Saved: part3_4cat_model_table_", filename_suffix, ".tex\n", sep = "")
  
  # Calculate SE inflation factors ---------------------------------------------
  cat("\n--- SE Inflation Analysis ---\n")
  
  inflation_3cat <- or_3cat_naive_95 %>%
    select(variable, naive_lower = conf_lower, naive_upper = conf_upper) %>%
    left_join(
      or_3cat_clusteredCY_95 %>%
        select(variable, clusterCY_lower = conf_lower, clusterCY_upper = conf_upper),
      by = "variable"
    ) %>%
    left_join(
      or_3cat_clusteredC_95 %>%
        select(variable, clusterC_lower = conf_lower, clusterC_upper = conf_upper),
      by = "variable"
    ) %>%
    mutate(
      ci_width_naive = naive_upper - naive_lower,
      ci_width_clusterCY = clusterCY_upper - clusterCY_lower,
      ci_width_clusterC = clusterC_upper - clusterC_lower,
      inflation_factor_CY = ci_width_clusterCY / ci_width_naive,
      inflation_factor_C = ci_width_clusterC / ci_width_naive
    )
  
  cat("3-category model SE inflation (City-Year):\n")
  cat("  Mean:", round(mean(inflation_3cat$inflation_factor_CY), 2), "\n")
  cat("  Median:", round(median(inflation_3cat$inflation_factor_CY), 2), "\n")
  cat("  Range:", round(min(inflation_3cat$inflation_factor_CY), 2), "-", 
      round(max(inflation_3cat$inflation_factor_CY), 2), "\n\n")
  
  cat("3-category model SE inflation (City):\n")
  cat("  Mean:", round(mean(inflation_3cat$inflation_factor_C), 2), "\n")
  cat("  Median:", round(median(inflation_3cat$inflation_factor_C), 2), "\n")
  cat("  Range:", round(min(inflation_3cat$inflation_factor_C), 2), "-", 
      round(max(inflation_3cat$inflation_factor_C), 2), "\n\n")
  
  # Create plots ---------------------------------------------------------------
  cat("--- Creating plots ---\n")
  
  # Helper function to create odds ratio plots
  create_or_plot <- function(or_data_95, or_data_90, output_prefix, se_label, subtitle_text, color_choice) {
    
    plot_data <- bind_rows(
      or_data_95 %>% mutate(ci_level = "95%"),
      or_data_90 %>% mutate(ci_level = "90%")
    ) %>%
      mutate(
        variable_clean = clean_var_names(variable),
        variable_clean = fct_reorder(variable_clean, or),
        ci_level = factor(ci_level, levels = c("90%", "95%"))
      )
    
    plot_data_95 <- plot_data %>% filter(ci_level == "95%")
    plot_data_90 <- plot_data %>% filter(ci_level == "90%")
    
    p <- ggplot() +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      # 90% CI - thicker linerange
      geom_linerange(
        data = plot_data_90,
        aes(y = variable_clean, xmin = conf_lower, xmax = conf_upper),
        linewidth = 1.5, color = color_choice
      ) +
      # 95% CI - thinner linerange
      geom_linerange(
        data = plot_data_95,
        aes(y = variable_clean, xmin = conf_lower, xmax = conf_upper),
        linewidth = 0.5, color = color_choice
      ) +
      # Points
      geom_point(
        data = plot_data_95,
        aes(x = or, y = variable_clean),
        size = 3, color = color_choice
      ) +
      labs(
        title = paste0("City-Level Determinants of State Capacity Signaling (", title_label, ")"),
        subtitle = subtitle_text,
        x = "Odds Ratio",
        y = NULL,
        caption = "Thick lines: 90% CI; Thin lines: 95% CI"
      ) +
      theme_minimal()
    
    filename <- paste0(output_prefix, "_", se_label, "_", filename_suffix, ".pdf")
    ggsave(file.path(output_dir, filename), 
           p, width = 10, height = 6, device = "pdf")
    cat("Saved:", filename, "\n")
  }
  
  # Plot 1: 3-category naive SEs
  create_or_plot(
    or_data_95 = or_3cat_naive_95,
    or_data_90 = or_3cat_naive_90,
    output_prefix = "part3_3cat",
    se_label = "naive",
    subtitle_text = "Administrative vs Coercive (Naive SEs)",
    color_choice = "steelblue"
  )
  
  # Plot 2: 3-category clustered SEs (City-Year)
  create_or_plot(
    or_data_95 = or_3cat_clusteredCY_95,
    or_data_90 = or_3cat_clusteredCY_90,
    output_prefix = "part3_3cat",
    se_label = "clusteredCY",
    subtitle_text = "Administrative vs Coercive (Cluster-Robust SEs by City-Year)",
    color_choice = "darkred"
  )
  
  # Plot 3: 3-category clustered SEs (City)
  create_or_plot(
    or_data_95 = or_3cat_clusteredC_95,
    or_data_90 = or_3cat_clusteredC_90,
    output_prefix = "part3_3cat",
    se_label = "clusteredC",
    subtitle_text = "Administrative vs Coercive (Cluster-Robust SEs by City)",
    color_choice = "darkgreen"
  )
  
  # Plot 4: 3-category comparison (all three)
  plot_data_3cat <- bind_rows(
    or_3cat_naive_95 %>% mutate(ci_level = "95%", se_type = "Naive SEs"),
    or_3cat_clusteredCY_95 %>% mutate(ci_level = "95%", se_type = "Clustered (City-Year)"),
    or_3cat_clusteredC_95 %>% mutate(ci_level = "95%", se_type = "Clustered (City)"),
    or_3cat_naive_90 %>% mutate(ci_level = "90%", se_type = "Naive SEs"),
    or_3cat_clusteredCY_90 %>% mutate(ci_level = "90%", se_type = "Clustered (City-Year)"),
    or_3cat_clusteredC_90 %>% mutate(ci_level = "90%", se_type = "Clustered (City)")
  ) %>%
    mutate(
      variable_clean = clean_var_names(variable),
      variable_clean = fct_reorder(variable_clean, or),
      se_type = factor(se_type, levels = c("Naive SEs", "Clustered (City-Year)", "Clustered (City)")),
      ci_level = factor(ci_level, levels = c("90%", "95%"))
    )
  
  plot_data_3cat_95 <- plot_data_3cat %>% filter(ci_level == "95%")
  plot_data_3cat_90 <- plot_data_3cat %>% filter(ci_level == "90%")
  
  plot_3cat_comparison <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    # 90% CI - thicker linerange
    geom_linerange(
      data = plot_data_3cat_90,
      aes(x = or, y = variable_clean, xmin = conf_lower, xmax = conf_upper,
          color = se_type),
      linewidth = 1.5,
      position = position_dodge(width = 0.6)
    ) +
    # 95% CI - thinner linerange
    geom_linerange(
      data = plot_data_3cat_95,
      aes(x = or, y = variable_clean, xmin = conf_lower, xmax = conf_upper,
          color = se_type),
      linewidth = 0.5,
      position = position_dodge(width = 0.6)
    ) +
    # Points
    geom_point(
      data = plot_data_3cat_95,
      aes(x = or, y = variable_clean, color = se_type, shape = se_type),
      size = 3,
      position = position_dodge(width = 0.6)
    ) +
    scale_color_manual(values = c("steelblue", "darkred", "darkgreen"), 
                       name = "Standard Error Type") +
    scale_shape_manual(values = c(16, 17, 15), name = "Standard Error Type") +
    labs(
      title = paste0("City-Level Determinants of State Capacity Signaling (", title_label, ")"),
      subtitle = "Administrative vs Coercive (Comparison)",
      x = "Odds Ratio",
      y = NULL,
      caption = "Thick lines: 90% CI; Thin lines: 95% CI"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("part3_3cat_comparison_", filename_suffix, ".pdf")), 
         plot_3cat_comparison, width = 12, height = 7, device = "pdf")
  cat("Saved: part3_3cat_comparison_", filename_suffix, ".pdf\n", sep = "")
  
  # Plots 5-8: 4-category Admin vs Coercive
  create_or_plot(
    or_data_95 = or_4cat_admin_naive_95,
    or_data_90 = or_4cat_admin_naive_90,
    output_prefix = "part3_4cat_admin",
    se_label = "naive",
    subtitle_text = "Administrative vs Coercive - 4-Category (Naive SEs)",
    color_choice = "steelblue"
  )
  
  create_or_plot(
    or_data_95 = or_4cat_admin_clusteredCY_95,
    or_data_90 = or_4cat_admin_clusteredCY_90,
    output_prefix = "part3_4cat_admin",
    se_label = "clusteredCY",
    subtitle_text = "Administrative vs Coercive - 4-Category (Cluster-Robust SEs by City-Year)",
    color_choice = "darkred"
  )
  
  create_or_plot(
    or_data_95 = or_4cat_admin_clusteredC_95,
    or_data_90 = or_4cat_admin_clusteredC_90,
    output_prefix = "part3_4cat_admin",
    se_label = "clusteredC",
    subtitle_text = "Administrative vs Coercive - 4-Category (Cluster-Robust SEs by City)",
    color_choice = "darkgreen"
  )
  
  # 4-category Admin comparison plot
  plot_data_4cat_admin <- bind_rows(
    or_4cat_admin_naive_95 %>% mutate(ci_level = "95%", se_type = "Naive SEs"),
    or_4cat_admin_clusteredCY_95 %>% mutate(ci_level = "95%", se_type = "Clustered (City-Year)"),
    or_4cat_admin_clusteredC_95 %>% mutate(ci_level = "95%", se_type = "Clustered (City)"),
    or_4cat_admin_naive_90 %>% mutate(ci_level = "90%", se_type = "Naive SEs"),
    or_4cat_admin_clusteredCY_90 %>% mutate(ci_level = "90%", se_type = "Clustered (City-Year)"),
    or_4cat_admin_clusteredC_90 %>% mutate(ci_level = "90%", se_type = "Clustered (City)")
  ) %>%
    mutate(
      variable_clean = clean_var_names(variable),
      variable_clean = fct_reorder(variable_clean, or),
      se_type = factor(se_type, levels = c("Naive SEs", "Clustered (City-Year)", "Clustered (City)")),
      ci_level = factor(ci_level, levels = c("90%", "95%"))
    )
  
  plot_data_4cat_admin_95 <- plot_data_4cat_admin %>% filter(ci_level == "95%")
  plot_data_4cat_admin_90 <- plot_data_4cat_admin %>% filter(ci_level == "90%")
  
  plot_4cat_admin_comparison <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_linerange(
      data = plot_data_4cat_admin_90,
      aes(x = or, y = variable_clean, xmin = conf_lower, xmax = conf_upper, color = se_type),
      linewidth = 1.5, position = position_dodge(width = 0.6)
    ) +
    geom_linerange(
      data = plot_data_4cat_admin_95,
      aes(x = or, y = variable_clean, xmin = conf_lower, xmax = conf_upper, color = se_type),
      linewidth = 0.5, position = position_dodge(width = 0.6)
    ) +
    geom_point(
      data = plot_data_4cat_admin_95,
      aes(x = or, y = variable_clean, color = se_type, shape = se_type),
      size = 3, position = position_dodge(width = 0.6)
    ) +
    scale_color_manual(values = c("steelblue", "darkred", "darkgreen"), name = "Standard Error Type") +
    scale_shape_manual(values = c(16, 17, 15), name = "Standard Error Type") +
    labs(
      title = paste0("City-Level Determinants - 4-Category (", title_label, ")"),
      subtitle = "Administrative vs Coercive (Comparison)",
      x = "Odds Ratio", y = NULL, caption = "Thick lines: 90% CI; Thin lines: 95% CI"
    ) +
    theme_minimal() + theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("part3_4cat_admin_comparison_", filename_suffix, ".pdf")), 
         plot_4cat_admin_comparison, width = 12, height = 7, device = "pdf")
  cat("Saved: part3_4cat_admin_comparison_", filename_suffix, ".pdf\n", sep = "")
  
  # Plots 9-12: 4-category Both vs Coercive
  create_or_plot(
    or_data_95 = or_4cat_both_naive_95,
    or_data_90 = or_4cat_both_naive_90,
    output_prefix = "part3_4cat_both",
    se_label = "naive",
    subtitle_text = "Coercive+Administrative vs Coercive - 4-Category (Naive SEs)",
    color_choice = "steelblue"
  )
  
  create_or_plot(
    or_data_95 = or_4cat_both_clusteredCY_95,
    or_data_90 = or_4cat_both_clusteredCY_90,
    output_prefix = "part3_4cat_both",
    se_label = "clusteredCY",
    subtitle_text = "Coercive+Administrative vs Coercive - 4-Category (Cluster-Robust SEs by City-Year)",
    color_choice = "darkred"
  )
  
  create_or_plot(
    or_data_95 = or_4cat_both_clusteredC_95,
    or_data_90 = or_4cat_both_clusteredC_90,
    output_prefix = "part3_4cat_both",
    se_label = "clusteredC",
    subtitle_text = "Coercive+Administrative vs Coercive - 4-Category (Cluster-Robust SEs by City)",
    color_choice = "darkgreen"
  )
  
  # 4-category Both comparison plot
  plot_data_4cat_both <- bind_rows(
    or_4cat_both_naive_95 %>% mutate(ci_level = "95%", se_type = "Naive SEs"),
    or_4cat_both_clusteredCY_95 %>% mutate(ci_level = "95%", se_type = "Clustered (City-Year)"),
    or_4cat_both_clusteredC_95 %>% mutate(ci_level = "95%", se_type = "Clustered (City)"),
    or_4cat_both_naive_90 %>% mutate(ci_level = "90%", se_type = "Naive SEs"),
    or_4cat_both_clusteredCY_90 %>% mutate(ci_level = "90%", se_type = "Clustered (City-Year)"),
    or_4cat_both_clusteredC_90 %>% mutate(ci_level = "90%", se_type = "Clustered (City)")
  ) %>%
    mutate(
      variable_clean = clean_var_names(variable),
      variable_clean = fct_reorder(variable_clean, or),
      se_type = factor(se_type, levels = c("Naive SEs", "Clustered (City-Year)", "Clustered (City)")),
      ci_level = factor(ci_level, levels = c("90%", "95%"))
    )
  
  plot_data_4cat_both_95 <- plot_data_4cat_both %>% filter(ci_level == "95%")
  plot_data_4cat_both_90 <- plot_data_4cat_both %>% filter(ci_level == "90%")
  
  plot_4cat_both_comparison <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_linerange(
      data = plot_data_4cat_both_90,
      aes(x = or, y = variable_clean, xmin = conf_lower, xmax = conf_upper, color = se_type),
      linewidth = 1.5, position = position_dodge(width = 0.6)
    ) +
    geom_linerange(
      data = plot_data_4cat_both_95,
      aes(x = or, y = variable_clean, xmin = conf_lower, xmax = conf_upper, color = se_type),
      linewidth = 0.5, position = position_dodge(width = 0.6)
    ) +
    geom_point(
      data = plot_data_4cat_both_95,
      aes(x = or, y = variable_clean, color = se_type, shape = se_type),
      size = 3, position = position_dodge(width = 0.6)
    ) +
    scale_color_manual(values = c("steelblue", "darkred", "darkgreen"), name = "Standard Error Type") +
    scale_shape_manual(values = c(16, 17, 15), name = "Standard Error Type") +
    labs(
      title = paste0("City-Level Determinants - 4-Category (", title_label, ")"),
      subtitle = "Coercive+Administrative vs Coercive (Comparison)",
      x = "Odds Ratio", y = NULL, caption = "Thick lines: 90% CI; Thin lines: 95% CI"
    ) +
    theme_minimal() + theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("part3_4cat_both_comparison_", filename_suffix, ".pdf")), 
         plot_4cat_both_comparison, width = 12, height = 7, device = "pdf")
  cat("Saved: part3_4cat_both_comparison_", filename_suffix, ".pdf\n", sep = "")
  
  # Save CSV files -------------------------------------------------------------
  or_3cat_naive_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_3cat_naive_95ci_", filename_suffix, ".csv")))
  or_3cat_clusteredCY_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_3cat_clusteredCY_95ci_", filename_suffix, ".csv")))
  or_3cat_clusteredC_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_3cat_clusteredC_95ci_", filename_suffix, ".csv")))
  or_3cat_naive_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_3cat_naive_90ci_", filename_suffix, ".csv")))
  or_3cat_clusteredCY_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_3cat_clusteredCY_90ci_", filename_suffix, ".csv")))
  or_3cat_clusteredC_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_3cat_clusteredC_90ci_", filename_suffix, ".csv")))
  inflation_3cat %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_3cat_se_inflation_", filename_suffix, ".csv")))
  
  or_4cat_admin_naive_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_admin_naive_95ci_", filename_suffix, ".csv")))
  or_4cat_admin_clusteredCY_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_admin_clusteredCY_95ci_", filename_suffix, ".csv")))
  or_4cat_admin_clusteredC_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_admin_clusteredC_95ci_", filename_suffix, ".csv")))
  or_4cat_admin_naive_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_admin_naive_90ci_", filename_suffix, ".csv")))
  or_4cat_admin_clusteredCY_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_admin_clusteredCY_90ci_", filename_suffix, ".csv")))
  or_4cat_admin_clusteredC_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_admin_clusteredC_90ci_", filename_suffix, ".csv")))
  
  or_4cat_both_naive_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_both_naive_95ci_", filename_suffix, ".csv")))
  or_4cat_both_clusteredCY_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_both_clusteredCY_95ci_", filename_suffix, ".csv")))
  or_4cat_both_clusteredC_95 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_both_clusteredC_95ci_", filename_suffix, ".csv")))
  or_4cat_both_naive_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_both_naive_90ci_", filename_suffix, ".csv")))
  or_4cat_both_clusteredCY_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_both_clusteredCY_90ci_", filename_suffix, ".csv")))
  or_4cat_both_clusteredC_90 %>% mutate(model = title_label) %>%
    write_csv(file.path(output_dir, paste0("part3_4cat_both_clusteredC_90ci_", filename_suffix, ".csv")))
  
  cat("\nModel fitting complete for:", title_label, "\n\n")
}

# ==============================================================================
# FIT MODELS WITH DIFFERENT SPECIFICATIONS
# ==============================================================================

predictor_labels_lookup <- tribble(
  ~predictor_name, ~clean_label,
  "city_lag1_grp_pc_yuan", "GDP per Capita",
  "city_lag1_population_total_10000", "Population",
  "city_lag1_distance_to_beijing_km", "Distance to Beijing",
  "city_lag1_urbanization_rate_pct", "Urbanization Rate",
  "city_lag1_education_avg_year", "Avg Yrs of Education",
  "city_lag1_public_exp_total_pc_10000yuan", "Public Expenditure per capita",
  "city_lag1_public_exp_publicservice_prop", "Public Service Expenditure %",
  "city_lag1_public_exp_pupblicsecurity_prop", "Public Security Expenditure %",
  "city_lag1_ethnic_minority_ratio", "Ethnic Minority Ratio",
  "city_lag1_unemployment_rate_urban", "Urban Unemployment Rate"
)

# (1) Development and State Capacity
fit_multinomial_models(
  predictors = c(
    "city_lag1_grp_pc_yuan",
    "city_lag1_population_total_10000",
    "city_lag1_distance_to_beijing_km",
    "city_lag1_urbanization_rate_pct",
    "city_lag1_education_avg_year"
  ),
  predictor_labels = predictor_labels_lookup,
  filename_suffix = "m1",
  title_label = "Development and State Capacity"
)

# (2) Public Expenditure
fit_multinomial_models(
  predictors = c(
    "city_lag1_public_exp_total_pc_10000yuan",
    "city_lag1_public_exp_publicservice_prop",
    "city_lag1_public_exp_pupblicsecurity_prop"
  ),
  filename_suffix = "m2",
  title_label = "Public Finance"
)


# (3) Social Stability
fit_multinomial_models(
  predictors = c(
    "city_lag1_ethnic_minority_ratio",
    "city_lag1_unemployment_rate_urban"
  ),
  filename_suffix = "m3",
  title_label = "Social Stability"
)


# (4) = (1) + (2)
fit_multinomial_models(
  predictors = c(
    "city_lag1_grp_pc_yuan",
    "city_lag1_population_total_10000",
    "city_lag1_distance_to_beijing_km",
    "city_lag1_urbanization_rate_pct",
    "city_lag1_education_avg_year",
    "city_lag1_public_exp_total_pc_10000yuan",
    "city_lag1_public_exp_publicservice_prop",
    "city_lag1_public_exp_pupblicsecurity_prop"
  ),
  filename_suffix = "m4",
  title_label = "Development, State Capacity, and Public Finance"
)


# (5) = (1) + (3)
fit_multinomial_models(
  predictors = c(
    "city_lag1_grp_pc_yuan",
    "city_lag1_population_total_10000",
    "city_lag1_distance_to_beijing_km",
    "city_lag1_urbanization_rate_pct",
    "city_lag1_education_avg_year",
    "city_lag1_ethnic_minority_ratio",
    "city_lag1_unemployment_rate_urban"
  ),
  filename_suffix = "m5",
  title_label = "Development, State Capacity, and Social Stability"
)


# (6) Full Model
fit_multinomial_models(
  predictors = c(
    "city_lag1_grp_pc_yuan",
    "city_lag1_population_total_10000",
    "city_lag1_distance_to_beijing_km",
    "city_lag1_urbanization_rate_pct",
    "city_lag1_education_avg_year",
    "city_lag1_public_exp_total_pc_10000yuan",
    "city_lag1_public_exp_publicservice_prop",
    "city_lag1_public_exp_pupblicsecurity_prop",
    "city_lag1_ethnic_minority_ratio",
    "city_lag1_unemployment_rate_urban"
  ),
  filename_suffix = "m6",
  title_label = "Full Model"
)



# ==============================================================================
# CREATE GROUPED FOREST PLOTS ACROSS ALL MODEL SPECIFICATIONS
# ==============================================================================

rm(list = ls())
output_dir <- "Figures_RR"

cat("\n===== Creating Grouped Forest Plots =====\n\n")

# Get all CSV files for each outcome/SE type combination
csv_files <- list.files(output_dir, pattern = "^part3_.+_[0-9]{2}ci_.+\\.csv$", full.names = TRUE)



# Helper function to clean variable names (remove backticks)
clean_var_names_plot <- function(vars) {
  str_remove_all(vars, "`")
}

# Function to create grouped forest plot
create_grouped_forest_plot <- function(outcome_type, se_type, outcome_comparison = NULL) {
  
  # For 4cat models, require outcome_comparison
  if (outcome_type %in% c("4cat_admin", "4cat_both") && is.null(outcome_comparison)) {
    stop("outcome_comparison must be specified for 4cat models")
  }
  
  # Read 95% CI files
  pattern_95 <- paste0("^part3_", outcome_type, "_", se_type, "_95ci_.+\\.csv$")
  files_95 <- list.files(output_dir, pattern = pattern_95, full.names = TRUE)
  
  # Read 90% CI files
  pattern_90 <- paste0("^part3_", outcome_type, "_", se_type, "_90ci_.+\\.csv$")
  files_90 <- list.files(output_dir, pattern = pattern_90, full.names = TRUE)
  
  if (length(files_95) == 0 || length(files_90) == 0) {
    cat("No files found for:", outcome_type, se_type, "\n")
    return(NULL)
  }
  
  cat("Creating grouped plot for:", outcome_type, se_type, "\n")
  cat("  95% CI files found:", length(files_95), "\n")
  cat("  90% CI files found:", length(files_90), "\n")
  
  # Read and combine all files
  combined_data_95 <- map_dfr(files_95, read_csv, col_types = cols(), show_col_types = FALSE)
  combined_data_90 <- map_dfr(files_90, read_csv, col_types = cols(), show_col_types = FALSE)
  
  # Clean variable names and prepare data
  plot_data_95 <- combined_data_95 %>%
    mutate(
      variable_clean = clean_var_names_plot(variable),
      model = factor(model, levels = unique(model)),
      ci_level = "95%"
    )
  
  plot_data_90 <- combined_data_90 %>%
    mutate(
      variable_clean = clean_var_names_plot(variable),
      model = factor(model, levels = unique(model)),
      ci_level = "90%"
    )
  
  # Create subtitle based on SE type and outcome
  if (se_type == "naive") {
    se_label <- "Naive SEs"
  } else if (se_type == "clusteredCY") {
    se_label <- "Cluster-Robust SEs (City-Year)"
  } else if (se_type == "clusteredC") {
    se_label <- "Cluster-Robust SEs (City)"
  }
  
  if (outcome_type == "3cat") {
    subtitle <- paste0("Administrative vs Coercive (", se_label, ")")
    outcome_label <- "3cat"
  } else if (outcome_type == "4cat_admin") {
    subtitle <- paste0(outcome_comparison, " (", se_label, ")")
    outcome_label <- "4cat_admin"
  } else if (outcome_type == "4cat_both") {
    subtitle <- paste0(outcome_comparison, " (", se_label, ")")
    outcome_label <- "4cat_both"
  }
  
  # Get number of unique models for color/shape scaling
  n_models <- length(unique(plot_data_95$model))
  
  # Create plot
  p <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    # 90% CI - thicker linerange
    geom_linerange(
      data = plot_data_90,
      aes(y = variable_clean, xmin = conf_lower, xmax = conf_upper, color = model),
      linewidth = 1,
      position = position_dodge(width = 0.4)
    ) +
    # 95% CI - thinner linerange
    geom_linerange(
      data = plot_data_95,
      aes(y = variable_clean, xmin = conf_lower, xmax = conf_upper, color = model),
      linewidth = 0.3,
      position = position_dodge(width = 0.4)
    ) +
    # Points
    geom_point(
      data = plot_data_95,
      aes(x = or, y = variable_clean, color = model, shape = model),
      size = 3,
      position = position_dodge(width = 0.4)
    ) +
    scale_color_manual(
      values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black")[1:n_models],
      name = "Model Specification"
    ) +
    scale_shape_manual(
      values = c(16, 17, 15, 18, 8, 4)[1:n_models],
      name = "Model Specification"
    ) +
    guides(
      color = guide_legend(ncol = 2, nrow = 3, byrow = TRUE),
      shape = guide_legend(ncol = 2, nrow = 3, byrow = TRUE)
    ) +
    labs(
      title = "City-Level Determinants of State Capacity Signaling",
      subtitle = subtitle,
      x = "Odds Ratio",
      y = NULL,
      caption = "Thick lines: 90% CI; Thin lines: 95% CI"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      # Larger axis text for readability
      axis.text.y = element_text(size = 14)
    )
  
  # Save plot with new naming convention
  filename <- paste0("part3_Z_", outcome_label, "_", se_type, "_grouped.pdf")
  ggsave(file.path(output_dir, filename), 
         p, width = 12, height = 8, device = "pdf")
  cat("Saved:", filename, "\n")
  
  return(p)
}

# Create all grouped forest plots
cat("\n===== Creating Grouped Forest Plots =====\n\n")

# 3-category models
create_grouped_forest_plot("3cat", "naive")
create_grouped_forest_plot("3cat", "clusteredCY")
create_grouped_forest_plot("3cat", "clusteredC")

# 4-category Admin vs Coercive models
create_grouped_forest_plot("4cat_admin", "naive", 
                           outcome_comparison = "Administrative vs Coercive")
create_grouped_forest_plot("4cat_admin", "clusteredCY", 
                           outcome_comparison = "Administrative vs Coercive")
create_grouped_forest_plot("4cat_admin", "clusteredC", 
                           outcome_comparison = "Administrative vs Coercive")

# 4-category Both vs Coercive models
create_grouped_forest_plot("4cat_both", "naive", 
                           outcome_comparison = "Coercive + Administrative vs Coercive")
create_grouped_forest_plot("4cat_both", "clusteredCY", 
                           outcome_comparison = "Coercive + Administrative vs Coercive")
create_grouped_forest_plot("4cat_both", "clusteredC", 
                           outcome_comparison = "Coercive + Administrative vs Coercive")

cat("\nPart III analysis complete.\n\n")

# ==============================================================================
# COMPLETION MESSAGE
# ==============================================================================

cat("\n===== Analysis Complete =====\n")
cat("All figures and tables saved to:", output_dir, "\n")
cat("\nSummary of outputs:\n")
cat("- Part 1: Prevalence plots\n")
cat("- Part 2: Event analysis\n")
cat("  * 3-category and 4-category models\n")
cat("  * Unified plots with naive and clustered SEs\n")
cat("  * LaTeX tables with corrected significance stars\n")
cat("  * 4-category: Admin vs Coercive and Both vs Coercive comparisons\n")
cat("- Part 3: City-level covariates\n")
cat("  * Multiple model specifications\n")
cat("  * 3-category and 4-category models\n")
cat("  * Odds ratio plots with naive and clustered SEs\n")
cat("  * LaTeX tables with corrected significance stars\n")
cat("  * 4-category: Admin vs Coercive and Both vs Coercive comparisons\n")
cat("  * Grouped forest plots comparing all model specifications\n")
