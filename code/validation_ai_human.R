################################################################################
# Purpose:
# This script compares AI-description ratings and human-description ratings.
#
# It creates four figures:
#
#   1. Videos where at least one AI rating is equal to or higher than
#      a human rating, by dimension.
#
#   2. Videos where at least one AI rating is higher than
#      a human rating, by dimension.
#
#   3. Videos where BOTH raters rated Agree or Strongly Agree,
#      with AI and Human bars side by side for each dimension.
#
#   4. Videos where AT LEAST ONE rater rated Agree or Strongly Agree,
#      with AI and Human bars side by side for each dimension.
################################################################################


################################################################################
# 1. Clean environment
################################################################################

rm(list = ls())
gc()
graphics.off()


################################################################################
# 2. Load required packages
################################################################################

# Install packages if needed:
# install.packages(c("tidyverse", "scales"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})


################################################################################
# 3. Define input file and output directory
################################################################################

data_path <- "data/description_validation_ai_human.rds"

# Output folder for plots and summary results
plot_dir <- "Figures_RR"

file_prefix <- "part6_"

output_file <- function(filename) {
  file.path(plot_dir, paste0(file_prefix, filename))
}

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}


################################################################################
# 4. Analysis settings
################################################################################

# Four evaluation dimensions
categories <- c(
  "coherence",
  "consistency",
  "fluency",
  "relevance"
)

category_labels <- c(
  "Coherence",
  "Consistency",
  "Fluency",
  "Relevance"
)

# Ratings considered positive
# +1 = Agree, +2 = Strongly Agree
positive_ratings <- c(1, 2)

# Allowed rating scale
allowed_ratings <- c(-2, -1, 0, 1, 2)

# Denominator for Figures 1 and 2:
# "at_least_one_each" = use video-dimension rows with at least one AI rating
#                       and at least one human rating.
# "complete_all_four" = require ai_eval1, ai_eval2, human_eval1, human_eval2
#                       all non-missing.
comparison_denominator <- "at_least_one_each"

# Denominator for Figures 3 and 4:
# "complete_ai_and_human_pairs" = require both AI raters and both human raters
#                                 for the same video-dimension.
# "source_valid_pairs"          = denominator is source-specific:
#                                 AI bar requires both AI raters;
#                                 Human bar requires both human raters.
positive_denominator <- "complete_ai_and_human_pairs"


################################################################################
# 5. Helper function: parse ratings
################################################################################

parse_likert_rating <- function(x) {
  
  x_chr <- as.character(x)
  x_chr <- str_replace_all(x_chr, "\u2212", "-") # Unicode minus to ASCII minus
  x_chr <- str_squish(x_chr)
  
  empty_like <- is.na(x_chr) |
    str_to_lower(x_chr) %in% c("", "na", "n/a", "null", "none", "missing")
  
  x_chr[empty_like] <- NA_character_
  x_low <- str_to_lower(x_chr)
  
  num <- suppressWarnings(parse_number(x_low))
  
  case_when(
    is.na(x_low) ~ NA_real_,
    str_detect(x_low, "strongly\\s+disagree") ~ -2,
    str_detect(x_low, "\\bdisagree\\b") ~ -1,
    str_detect(x_low, "neutral|neither") ~ 0,
    str_detect(x_low, "strongly\\s+agree") ~ 2,
    str_detect(x_low, "\\bagree\\b") ~ 1,
    !is.na(num) ~ as.numeric(num),
    TRUE ~ NA_real_
  )
}


################################################################################
# 6. Load dataset
################################################################################

validation_df <- readRDS(data_path)

# Remove accidental leading/trailing spaces in column names
names(validation_df) <- str_squish(names(validation_df))

# Add stable row id
validation_df <- validation_df %>%
  mutate(row_id_temp = row_number())

validation_df <- validation_df %>%
  mutate(video_id_temp = paste0("row_", row_id_temp))
stopifnot("human_available" %in% names(validation_df))


################################################################################
# 7. Define required rating columns
################################################################################

ai_eval1_cols <- paste0("ai_eval1_", categories)
ai_eval2_cols <- paste0("ai_eval2_", categories)

human_eval1_cols <- paste0("human_eval1_", categories)
human_eval2_cols <- paste0("human_eval2_", categories)

required_cols <- c(
  ai_eval1_cols,
  ai_eval2_cols,
  human_eval1_cols,
  human_eval2_cols
)

missing_cols <- setdiff(required_cols, names(validation_df))

if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from the dataset:\n",
    paste(missing_cols, collapse = ", ")
  )
}

if (!comparison_denominator %in% c("at_least_one_each", "complete_all_four")) {
  stop("comparison_denominator must be 'at_least_one_each' or 'complete_all_four'.")
}

if (!positive_denominator %in% c("complete_ai_and_human_pairs", "source_valid_pairs")) {
  stop("positive_denominator must be 'complete_ai_and_human_pairs' or 'source_valid_pairs'.")
}


################################################################################
# 8. Reshape ratings to long format
################################################################################

ratings_long <- validation_df %>%
  select(
    row_id_temp,
    video_id_temp,
    human_available,
    all_of(required_cols)
  ) %>%
  pivot_longer(
    cols = all_of(required_cols),
    names_to = c("source", "rater_num", "category"),
    names_pattern = "^(ai|human)_eval([12])_(.*)$",
    values_to = "rating_raw"
  ) %>%
  mutate(
    rater = paste0("eval", rater_num),
    rating_num = parse_likert_rating(rating_raw)
  )

# Warn about ratings outside the allowed scale
invalid_ratings <- ratings_long %>%
  filter(
    !is.na(rating_num),
    !(rating_num %in% allowed_ratings)
  )

if (nrow(invalid_ratings) > 0) {
  warning(
    "Some parsed ratings are outside the allowed -2 to +2 scale. ",
    "They will be treated as missing."
  )
}

ratings_long <- ratings_long %>%
  mutate(
    rating_num = if_else(
      rating_num %in% allowed_ratings,
      rating_num,
      NA_real_
    )
  )


################################################################################
# 9. Reshape ratings to one row per video x dimension
################################################################################

ratings_wide <- ratings_long %>%
  select(
    row_id_temp,
    video_id_temp,
    human_available,
    category,
    source,
    rater,
    rating_num
  ) %>%
  pivot_wider(
    id_cols = c(
      row_id_temp,
      video_id_temp,
      human_available,
      category
    ),
    names_from = c(source, rater),
    values_from = rating_num,
    names_glue = "{source}_{rater}"
  ) %>%
  mutate(
    category = factor(
      category,
      levels = categories,
      labels = category_labels
    )
  )


################################################################################
# 10. Create comparison variables
################################################################################

comparison_wide <- ratings_wide %>%
  rowwise() %>%
  mutate(
    n_ai_ratings = sum(!is.na(c_across(c(ai_eval1, ai_eval2)))),
    n_human_ratings = sum(!is.na(c_across(c(human_eval1, human_eval2)))),
    
    ai_max = {
      vals <- c_across(c(ai_eval1, ai_eval2))
      if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
    },
    
    human_min = {
      vals <- c_across(c(human_eval1, human_eval2))
      if (all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE)
    }
  ) %>%
  ungroup() %>%
  mutate(
    in_human_subset = human_available | n_human_ratings > 0,
    
    valid_comparison_at_least_one_each =
      n_ai_ratings >= 1 & n_human_ratings >= 1,
    
    valid_comparison_complete_all_four =
      n_ai_ratings == 2 & n_human_ratings == 2,
    
    valid_comparison =
      if (comparison_denominator == "at_least_one_each") {
        valid_comparison_at_least_one_each
      } else {
        valid_comparison_complete_all_four
      },
    
    # Figure 1 condition:
    # There exists at least one AI rating and one human rating such that AI >= Human.
    # Equivalent to max(AI ratings) >= min(Human ratings).
    ai_ge_any_human = ai_max >= human_min,
    
    # Figure 2 condition:
    # There exists at least one AI rating and one human rating such that AI > Human.
    # Equivalent to max(AI ratings) > min(Human ratings).
    ai_gt_any_human = ai_max > human_min,
    
    valid_ai_pair = n_ai_ratings == 2,
    valid_human_pair = n_human_ratings == 2
  )


################################################################################
# 11. Summary for Figures 1 and 2
################################################################################

comparison_analysis <- comparison_wide %>%
  filter(
    in_human_subset,
    valid_comparison
  )

comparison_denominator_note <- if (comparison_denominator == "at_least_one_each") {
  "Human-description videos with at least one AI rating and at least one human rating"
} else {
  "Human-description videos with both AI ratings and both human ratings"
}

comparison_summary <- comparison_analysis %>%
  group_by(category, .drop = FALSE) %>%
  summarise(
    denominator_n = n(),
    
    ai_ge_human_n = sum(ai_ge_any_human, na.rm = TRUE),
    ai_gt_human_n = sum(ai_gt_any_human, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    ai_ge_human_pct = if_else(
      denominator_n > 0,
      100 * ai_ge_human_n / denominator_n,
      NA_real_
    ),
    
    ai_gt_human_pct = if_else(
      denominator_n > 0,
      100 * ai_gt_human_n / denominator_n,
      NA_real_
    ),
    
    denominator_rule = comparison_denominator_note,
    comparison_rule = "Any-pair rule: max(AI ratings) compared with min(Human ratings)"
  )

write_csv(
  comparison_summary,
  output_file("ai_human_rating_comparison_summary.csv")
)


################################################################################
# 12. Summary for Figures 3 and 4
################################################################################

positive_base <- comparison_wide %>%
  filter(in_human_subset)

ai_positive <- positive_base %>%
  transmute(
    row_id_temp,
    video_id_temp,
    category,
    source = "AI",
    
    source_valid_pair = valid_ai_pair,
    both_sources_complete_pairs = valid_ai_pair & valid_human_pair,
    
    both_raters_positive =
      ai_eval1 %in% positive_ratings &
      ai_eval2 %in% positive_ratings,
    
    at_least_one_rater_positive =
      ai_eval1 %in% positive_ratings |
      ai_eval2 %in% positive_ratings
  )

human_positive <- positive_base %>%
  transmute(
    row_id_temp,
    video_id_temp,
    category,
    source = "Human",
    
    source_valid_pair = valid_human_pair,
    both_sources_complete_pairs = valid_ai_pair & valid_human_pair,
    
    both_raters_positive =
      human_eval1 %in% positive_ratings &
      human_eval2 %in% positive_ratings,
    
    at_least_one_rater_positive =
      human_eval1 %in% positive_ratings |
      human_eval2 %in% positive_ratings
  )

positive_long <- bind_rows(ai_positive, human_positive) %>%
  mutate(
    source = factor(source, levels = c("AI", "Human"))
  )

if (positive_denominator == "complete_ai_and_human_pairs") {
  
  positive_analysis <- positive_long %>%
    filter(both_sources_complete_pairs)
  
  positive_denominator_note <-
    "Same video-dimension rows with both AI raters and both human raters"
  
} else {
  
  positive_analysis <- positive_long %>%
    filter(source_valid_pair)
  
  positive_denominator_note <-
    "Source-specific valid pairs: AI requires both AI raters; Human requires both human raters"
}

positive_summary <- positive_analysis %>%
  group_by(category, source, .drop = FALSE) %>%
  summarise(
    denominator_n = n(),
    
    both_raters_n = sum(
      both_raters_positive,
      na.rm = TRUE
    ),
    
    at_least_one_rater_n = sum(
      at_least_one_rater_positive,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    both_raters_pct = if_else(
      denominator_n > 0,
      100 * both_raters_n / denominator_n,
      NA_real_
    ),
    
    at_least_one_rater_pct = if_else(
      denominator_n > 0,
      100 * at_least_one_rater_n / denominator_n,
      NA_real_
    ),
    
    denominator_rule = positive_denominator_note
  )

write_csv(
  positive_summary,
  output_file("ai_vs_human_positive_rating_summary.csv")
)


################################################################################
# 13. Plotting functions
################################################################################

ai_plot_label <- "LVLM"

plot_percentage_bar <- function(plot_df, plot_title, fill_color) {
  
  plot_df <- plot_df %>%
    mutate(
      category = factor(as.character(category), levels = category_labels),
      pct_plot = replace_na(pct, 0),
      label = if_else(
        denominator_n > 0,
        paste0(
          number(pct, accuracy = 0.1),
          "%\n(",
          numerator_n,
          "/",
          denominator_n,
          ")"
        ),
        "No valid\nratings"
      )
    )
  
  ggplot(plot_df, aes(x = category, y = pct_plot)) +
    geom_col(
      width = 0.65,
      fill = fill_color
    ) +
    geom_text(
      aes(label = label),
      vjust = -0.35,
      size = 4,
      lineheight = 0.9
    ) +
    scale_y_continuous(
      limits = c(0, 110),
      breaks = seq(0, 100, by = 20),
      labels = label_percent(scale = 1)
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Percentage of videos"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    )
}


plot_grouped_percentage_bar <- function(plot_df, plot_title) {
  
  dodge_width <- 0.75
  
  plot_df <- plot_df %>%
    mutate(
      category = factor(as.character(category), levels = category_labels),
      source = factor(
        as.character(source),
        levels = c("AI", "Human"),
        labels = c(ai_plot_label, "Human")
      ),
      pct_plot = replace_na(pct, 0),
      label = if_else(
        denominator_n > 0,
        paste0(
          number(pct, accuracy = 0.1),
          "%\n(",
          numerator_n,
          "/",
          denominator_n,
          ")"
        ),
        "No valid\nratings"
      )
    )
  
  ggplot(plot_df, aes(x = category, y = pct_plot, fill = source)) +
    geom_col(
      position = position_dodge(width = dodge_width),
      width = 0.68,
      color = "white"
    ) +
    geom_text(
      aes(label = label),
      position = position_dodge(width = dodge_width),
      vjust = -0.30,
      size = 3.5,
      lineheight = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "LVLM" = "#2C7FB8",
        "Human" = "#F28E2B"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, 110),
      breaks = seq(0, 100, by = 20),
      labels = label_percent(scale = 1)
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Percentage of videos",
      fill = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    )
}


################################################################################
# 14. Figure 1:
# Videos where at least one AI rating is equal to or higher than a human rating
################################################################################

plot_ai_ge_human_df <- comparison_summary %>%
  transmute(
    category,
    denominator_n,
    numerator_n = ai_ge_human_n,
    pct = ai_ge_human_pct
  )

plot_ai_ge_human <- plot_percentage_bar(
  plot_df = plot_ai_ge_human_df,
  plot_title = "Videos Where At Least One LVLM Rating Is Equal to or Higher Than a Human Rating",
  fill_color = "#756BB1"
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 10.5
    )
  )

print(plot_ai_ge_human)


################################################################################
# 15. Figure 2:
# Videos where at least one AI rating is higher than a human rating
################################################################################

plot_ai_gt_human_df <- comparison_summary %>%
  transmute(
    category,
    denominator_n,
    numerator_n = ai_gt_human_n,
    pct = ai_gt_human_pct
  )

plot_ai_gt_human <- plot_percentage_bar(
  plot_df = plot_ai_gt_human_df,
  plot_title = "Videos Where At Least One LVLM Rating Is Higher Than a Human Rating",
  fill_color = "#E6550D"
)

print(plot_ai_gt_human)


################################################################################
# 16. Figure 3:
# Videos where at least one rater rated Agree or Strongly Agree
################################################################################

plot_at_least_one_positive_df <- positive_summary %>%
  transmute(
    category,
    source,
    denominator_n,
    numerator_n = at_least_one_rater_n,
    pct = at_least_one_rater_pct
  )

plot_at_least_one_positive <- plot_grouped_percentage_bar(
  plot_df = plot_at_least_one_positive_df,
  plot_title = "Videos Where At Least One Rater Rated Agree or Strongly Agree"
)

print(plot_at_least_one_positive)


################################################################################
# 17. Save plots
################################################################################

save_plot_png_pdf <- function(plot_object, filename_stub, width = 8, height = 5.5) {
  
  ggsave(
    filename = output_file(paste0(filename_stub, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = output_file(paste0(filename_stub, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    bg = "white"
  )
}

save_plot_png_pdf(
  plot_ai_ge_human,
  "fig1_pct_ai_rating_equal_or_higher_than_human_rating"
)

save_plot_png_pdf(
  plot_ai_gt_human,
  "fig2_pct_ai_rating_higher_than_human_rating"
)

save_plot_png_pdf(
  plot_at_least_one_positive,
  "fig4_pct_at_least_one_rater_agree_or_strongly_agree_ai_vs_human"
)


################################################################################
# 19. Save one combined 2 x 2 figure if patchwork is installed
################################################################################

if (requireNamespace("patchwork", quietly = TRUE)) {
  
  combined_plot <- patchwork::wrap_plots(
    plot_ai_ge_human,
    plot_ai_gt_human,
    plot_at_least_one_positive,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = "LVLM and Human Description Rating Comparisons"
    )
  
  save_plot_png_pdf(
    combined_plot,
    "combined_ai_human_rating_comparison_figures",
    width = 16,
    height = 11
  )
}


################################################################################
# End of script
################################################################################
