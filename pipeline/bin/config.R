#' Load environment variables from .env file
#' @param env_file Path to .env file
load_dot_env <- function(env_file = ".env") {
  if (file.exists(env_file)) {
    lines <- readLines(env_file, warn = FALSE)
    lines <- lines[!grepl("^\\s*(#|$)", lines)] # Remove comments and empty lines
    for (line in lines) {
      parts <- strsplit(line, "=")[[1]]
      if (length(parts) >= 2) {
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = "="))
        # Remove quotes if present
        value <- gsub("^['\"]|['\"]$", "", value)
        env_list <- list(value)
        names(env_list) <- key
        do.call(Sys.setenv, env_list)
      }
    }
  }
}

# Try to load .env from project root
# 1. Check GWAS_PROJECT_ROOT if already set
# 2. Check current directory
# 3. Check one level up (if run from bin/ or src/)
# 4. Check two levels up (if run from subfolders)
project_root <- Sys.getenv("GWAS_PROJECT_ROOT")

search_paths <- c(
    if (project_root != "") project_root else NULL,
    ".",
    "..",
    "../.."
)

for (p in search_paths) {
    env_p <- file.path(p, ".env")
    if (file.exists(env_p)) {
        load_dot_env(env_p)
        break
    }
}

get_env <- function(key, default = NULL) {
  val <- Sys.getenv(key)
  if (val == "") return(default)
  return(val)
}
