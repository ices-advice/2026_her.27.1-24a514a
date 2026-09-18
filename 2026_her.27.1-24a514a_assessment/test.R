#------------------------------------------------------------------------------#
# test.R
#
# This is a quality assurance script!
# It tests if files are reproducible on github Workflow
#
# It is possible to rename the folders named from XXXX to XXXX_reference, and 
# run this script. This can be used when optimising 
#
# TODO: check if it possible to test if Rda are identical
#------------------------------------------------------------------------------#



        
reference_suffix <- "_reference"
check_dirs       <- c("data", "model", "output", "report","forecast")

# Files to skip entirely (path relative to the check_dirs entry, e.g. "forecast/...")
ignore_files <- c(
  "forecast/tables/tab4.8.2.1b_alternative.txt"
)

all_failures <- list()

compare_dir <- function(dir_name) {
  reference_dir <- paste0(dir_name, reference_suffix)
  new_dir       <- dir_name
  failures      <- character(0)

  if (!dir.exists(reference_dir)) {
    cat(sprintf("No reference found for %s/ — skipping.\n", dir_name))
    return(failures)
  }

  cat(sprintf("\nChecking %s/ ...\n", dir_name))

  # Compare CSV files
  csv_files <- list.files(reference_dir, pattern = "\\.csv$", full.names = F,recursive = T)
  csv_files <- csv_files[!file.path(dir_name, csv_files) %in% ignore_files]
  for (f in csv_files) {
    ref    <- read.csv(file.path(reference_dir, f), check.names = FALSE)
    new    <- read.csv(file.path(new_dir, f),       check.names = FALSE)
    result <- all.equal(ref, new)
    if (!isTRUE(result)) {
      failures <- c(failures, paste0(dir_name, "/", f, ": ", paste(result, collapse = "; ")))
    } else {
      cat(sprintf("  ✓ %s\n", f))
    }
  }

  # Compare TXT files
  txt_files <- list.files(reference_dir, pattern = "\\.txt$", full.names = FALSE,recursive = T)
  txt_files <- txt_files[!file.path(dir_name, txt_files) %in% ignore_files]
  for (f in txt_files) {
    ref <- readLines(file.path(reference_dir, f))
    new <- readLines(file.path(new_dir,       f))
    if (!identical(ref, new)) {
      failures <- c(failures, paste0(dir_name, "/", f, ": content differs"))
    } else {
      cat(sprintf("  ✓ %s\n", f))
    }
  }

  failures
}

# Run checks for all directories
for (d in check_dirs) {
  all_failures[[d]] <- compare_dir(d)
}

# Flatten failures list
failures_flat <- unlist(all_failures)

# Final report — crash here so all failures are shown before stopping
cat("\n========================================\n")
if (length(failures_flat) == 0) {
  cat("All results are reproducible!\n")
} else {
  cat(sprintf("Reproducibility check FAILED (%d issue(s)):\n", length(failures_flat)))
  cat(paste0("  ✗ ", failures_flat, "\n"))
  cat("========================================\n")
  stop("Reproducibility check failed — see details above.")
}
cat("========================================\n")