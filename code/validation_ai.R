################################################################################
# Purpose:
# This script compares validation results from two raters/evaluators.
# It calculates and plots:
#   1. Percentage of videos where BOTH raters rated Agree (1) or Strongly Agree (2)
#      across four evaluation categories:
#        - coherence
#        - consistency
#        - fluency
#        - relevance
#
#   2. Percentage of videos where AT LEAST ONE rater rated Agree (1) or
#      Strongly Agree (2) across the same four categories.
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
# 3. Define input file and plot output directory
################################################################################

# Update this path to your dataset
data_path <- "data/description_validation_ai.rds"

# Output folder for plots and summary results
plot_dir <- "Figures_RR"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}


################################################################################
# 4. Load dataset
################################################################################

validation_df <- readRDS(data_path)

# Remove accidental leading/trailing spaces in column names
names(validation_df) <- trimws(names(validation_df))


################################################################################
# 5. Define columns and analysis settings
################################################################################

# Four validation categories
categories <- c(
  "coherence",
  "consistency",
  "fluency",
  "relevance"
)

# Evaluator columns
eval1_cols <- paste0("eval1_", categories)
eval2_cols <- paste0("eval2_", categories)

required_cols <- c(eval1_cols, eval2_cols)

# Check that all required columns exist
missing_cols <- setdiff(required_cols, names(validation_df))

if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from the dataset:\n",
    paste(missing_cols, collapse = ", ")
  )
}

# Ratings considered positive
# 1 = Agree, 2 = Strongly Agree
positive_ratings <- c(1, 2)

# Denominator option:
# "valid_pairs" = use only videos where both raters have non-missing ratings
# "all_videos"  = use all rows/videos in the dataset
denominator <- "valid_pairs"


################################################################################
# 6. Prepare data for comparison
################################################################################

validation_long <- validation_df %>%
  mutate(video_id_temp = row_number()) %>%
  select(video_id_temp, all_of(required_cols)) %>%
  pivot_longer(
    cols = all_of(required_cols),
    names_to = c("rater", "category"),
    names_pattern = "^(eval[12])_(.*)$",
    values_to = "rating_raw"
  ) %>%
  mutate(
    category = factor(
      category,
      levels = categories,
      labels = c("Coherence", "Consistency", "Fluency", "Relevance")
    ),
    
    # Convert ratings to numeric.
    # This works for values like 1, 2, "1", "2", "Agree (1)", etc.
    rating_num = suppressWarnings(parse_number(as.character(rating_raw))),
    
    # Also handle pure text values if present
    rating_text = str_to_lower(str_squish(as.character(rating_raw))),
    rating_num = case_when(
      is.na(rating_num) & rating_text == "agree" ~ 1,
      is.na(rating_num) & rating_text == "strongly agree" ~ 2,
      TRUE ~ rating_num
    ),
    
    agree_or_strongly_agree = rating_num %in% positive_ratings
  )


################################################################################
# 7. Create paired evaluator data
################################################################################

paired_validation <- validation_long %>%
  select(
    video_id_temp,
    rater,
    category,
    rating_num,
    agree_or_strongly_agree
  ) %>%
  pivot_wider(
    names_from = rater,
    values_from = c(rating_num, agree_or_strongly_agree),
    names_glue = "{.value}_{rater}"
  ) %>%
  mutate(
    valid_pair = !is.na(rating_num_eval1) & !is.na(rating_num_eval2),
    
    agree_or_strongly_agree_eval1 =
      coalesce(agree_or_strongly_agree_eval1, FALSE),
    
    agree_or_strongly_agree_eval2 =
      coalesce(agree_or_strongly_agree_eval2, FALSE),
    
    both_raters_agree_or_strongly_agree =
      agree_or_strongly_agree_eval1 & agree_or_strongly_agree_eval2,
    
    at_least_one_rater_agree_or_strongly_agree =
      agree_or_strongly_agree_eval1 | agree_or_strongly_agree_eval2
  )


################################################################################
# 8. Apply denominator rule
################################################################################

if (denominator == "valid_pairs") {
  
  analysis_df <- paired_validation %>%
    filter(valid_pair)
  
  denominator_note <- "Videos with non-missing ratings from both raters"
  
} else if (denominator == "all_videos") {
  
  analysis_df <- paired_validation
  
  denominator_note <- "All videos in dataset; missing ratings treated as not positive"
  
} else {
  
  stop("denominator must be either 'valid_pairs' or 'all_videos'")
}


################################################################################
# 9. Calculate summary percentages
################################################################################

summary_results <- analysis_df %>%
  group_by(category, .drop = FALSE) %>%
  summarise(
    denominator_n = n(),
    
    both_raters_n = sum(
      both_raters_agree_or_strongly_agree,
      na.rm = TRUE
    ),
    
    at_least_one_rater_n = sum(
      at_least_one_rater_agree_or_strongly_agree,
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
    )
  )

# Save summary table
write_csv(
  summary_results,
  file.path(plot_dir, "validation_agreement_summary.csv")
)


################################################################################
# 10. Prepare plotting function
################################################################################

plot_percentage_bar <- function(plot_df, plot_title, fill_color) {
  
  plot_df <- plot_df %>%
    mutate(
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
      size = 4
    ) +
    scale_y_continuous(
      limits = c(0, 105),
      breaks = seq(0, 100, by = 20),
      labels = label_percent(scale = 1)
    ) +
    labs(
      title = plot_title,
      # subtitle = denominator_note,
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


################################################################################
# 11. Plot:
# Percentage of videos where AT LEAST ONE rater rated Agree or Strongly Agree
################################################################################

plot_at_least_one_df <- summary_results %>%
  transmute(
    category,
    denominator_n,
    numerator_n = at_least_one_rater_n,
    pct = at_least_one_rater_pct
  )

plot_at_least_one <- plot_percentage_bar(
  plot_df = plot_at_least_one_df,
  plot_title = "Videos Where At Least One Rater Rated Agree or Strongly Agree",
  fill_color = "#41AB5D"
)

print(plot_at_least_one)


################################################################################
# 12. Save plots
################################################################################


ggsave(
  filename = file.path(
    plot_dir,
    "plot2_pct_at_least_one_rater_agree_or_strongly_agree.png"
  ),
  plot = plot_at_least_one,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)


################################################################################
# 14. Save plots as PDF as well
################################################################################

ggsave(
  filename = file.path(
    plot_dir,
    "plot2_pct_at_least_one_rater_agree_or_strongly_agree.pdf"
  ),
  plot = plot_at_least_one,
  width = 8,
  height = 5.5,
  bg = "white"
)


################################################################################
# End of script
################################################################################
