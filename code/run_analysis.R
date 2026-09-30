run_analysis <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  stopifnot(file.exists(file.path(root, "data", "analysis.rds")),
            file.exists(file.path(root, "code", "analysis_events_city.R")))
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") {
    "Rscript.exe"
  } else {
    "Rscript"
  })
  stopifnot(file.exists(rscript))

  work <- tempfile("douyin_replication_")
  dir.create(work)
  dir.create(file.path(work, "data"))
  file.copy(list.files(file.path(root, "data"), full.names = TRUE),
            file.path(work, "data"))
  out <- file.path(root, "Figures_Tables")
  dir.create(out, showWarnings = FALSE)

  run_script <- function(filename) {
    log <- file.path(work, paste0(filename, ".log"))
    old_wd <- getwd()
    old_lc_all <- Sys.getenv("LC_ALL", unset = "")
    old_lc_ctype <- Sys.getenv("LC_CTYPE", unset = "")
    on.exit({
      setwd(old_wd)
      Sys.setenv(LC_ALL = old_lc_all, LC_CTYPE = old_lc_ctype)
    }, add = TRUE)
    setwd(work)
    if (.Platform$OS.type == "windows") {
      Sys.setenv(LC_ALL = ".UTF-8", LC_CTYPE = ".UTF-8")
    }
    status <- system2(
      rscript,
      c("--vanilla", shQuote(file.path(root, "code", filename))),
      stdout = log, stderr = log
    )
    if (!identical(status, 0L)) {
      stop(filename, " failed:\n", paste(tail(readLines(log, warn = FALSE), 25L),
                                          collapse = "\n"))
    }
    invisible(NULL)
  }

  files <- list(
    analysis_events_city.R = c(
      "part1_prevalence_by_time_3cat.pdf" = "Figure_3.pdf",
      "part1_prevalence_by_city_3cat.pdf" = "Figure_4.pdf",
      "part2_3cat_unified_clustered.pdf" = "Figure_5.pdf",
      "part3_Z_3cat_clusteredC_grouped.pdf" = "Figure_6.pdf",
      "part1_prevalence_by_time_4cat.pdf" = "Appendix_Figure_4.pdf",
      "part1_prevalence_by_city_4cat.pdf" = "Appendix_Figure_5.pdf",
      "part2_3cat_models_clustered_se.tex" = "Table_A3.tex",
      "part3_3cat_model_table_m1.tex" = "Table_A5.tex",
      "part3_3cat_model_table_m2.tex" = "Table_A6.tex",
      "part3_3cat_model_table_m3.tex" = "Table_A7.tex",
      "part3_3cat_model_table_m4.tex" = "Table_A8.tex",
      "part3_3cat_model_table_m5.tex" = "Table_A9.tex",
      "part3_3cat_model_table_m6.tex" = "Table_A10.tex",
      "part3_4cat_model_table_m1.tex" = "Table_A13.tex",
      "part3_4cat_model_table_m2.tex" = "Table_A14.tex",
      "part3_4cat_model_table_m3.tex" = "Table_A15.tex",
      "part3_4cat_model_table_m4.tex" = "Table_A16.tex",
      "part3_4cat_model_table_m6.tex" = "Table_A17.tex",
      "part2_4cat_models_clustered_se.tex" = "Table_A18.tex"
    ),
    analysis_engagement.R = c(
      "figure_nb_pred_engagement_admin_vs_coercive.pdf" = "Figure_7.pdf",
      "table_nb_engagement_admin_vs_coercive.tex" = "Table_A11.tex",
      "table_nb_engagement_admin_vs_coercive_with_followers.tex" = "Table_A12.tex",
      "table_nb_engagement_4cat_no_followers.tex" = "Table_A19.tex",
      "table_nb_engagement_4cat_with_followers.tex" = "Table_A20.tex"
    ),
    validation_ai.R = c(
      "plot2_pct_at_least_one_rater_agree_or_strongly_agree.pdf" = "Appendix_Figure_2.pdf"
    ),
    validation_ai_human.R = c(
      "part6_fig4_pct_at_least_one_rater_agree_or_strongly_agree_ai_vs_human.pdf" =
        "Appendix_Figure_3.pdf"
    )
  )

  for (script in names(files)) {
    run_script(script)
    for (source_name in names(files[[script]])) {
      source_path <- file.path(work, "Figures_RR", source_name)
      target_path <- file.path(out, unname(files[[script]][[source_name]]))
      if (!file.exists(source_path)) {
        stop("Expected output missing from ", script, ": ", source_name)
      }
      if (!file.copy(source_path, target_path, overwrite = TRUE)) {
        stop("Could not copy ", source_name, " to Figures_Tables")
      }
    }
  }
  invisible(unlist(files, use.names = FALSE))
}
