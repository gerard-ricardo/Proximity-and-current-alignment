
#How to populate this project
#1) Create folder called 'R_scripts'. Copy all script from Palau project to this folder manually.
#2) Run below code
#3) Delete 'R_scripts'
#3) Push to Github






library(fs)

input_folder <- "C:/Users/gerar/Documents/1_R_projects_local/Proximity-and-current-alignment/R_scripts"
output_folder <- "C:/Users/gerar/Documents/1_R_projects_local/Proximity-and-current-alignment/R_scripts3"

remove_single_hash_comment <- function(line) {
  # Find positions of single and double quotes
  double_quotes <- gregexpr('"', line)[[1]]
  single_quotes <- gregexpr("'", line)[[1]]
  
  # Combine and sort quote positions with type labels
  quote_positions <- c(double_quotes, single_quotes)
  if (all(quote_positions == -1)) return(sub("\\s*#.*$", "", line)) # no quotes, remove comment
  quote_positions <- sort(quote_positions[quote_positions != -1])
  
  # Get corresponding quote types for each position
  quote_types <- c(rep('"', length(double_quotes[double_quotes != -1])),
                   rep("'", length(single_quotes[single_quotes != -1])))
  quote_types <- quote_types[order(c(double_quotes[double_quotes != -1], single_quotes[single_quotes != -1]))]
  
  hash_pos <- regexpr("#", line)[1]
  if (hash_pos == -1) return(line) # no hash, keep line
  
  inside_quotes <- FALSE
  # Iterate over quotes in pairs
  for (i in seq(1, length(quote_positions), 2)) {
    start <- quote_positions[i]
    end <- ifelse(i + 1 <= length(quote_positions), quote_positions[i + 1], nchar(line) + 1)
    if (hash_pos > start && hash_pos < end) {
      inside_quotes <- TRUE
      break
    }
  }
  if (inside_quotes) return(line) # keep line if # inside quotes
  sub("\\s*#.*$", "", line) # else remove comment
}


clean_project_code <- function(input_dir, output_dir) {
  dir_create(output_dir)
  r_files <- dir_ls(input_dir, recurse = TRUE, glob = "*.R")
  for (file in r_files) {
    rel_path <- path_rel(file, start = input_dir)
    out_file <- path(output_dir, rel_path)
    dir_create(path_dir(out_file))
    lines <- readLines(file)
    cleaned <- character(length(lines))
    j <- 1
    for (i in seq_along(lines)) {
      line <- lines[i]
      if (grepl("^\\s*##", line) || grepl("^\\s*#.*-{4,}", line)) {
        if (j > 1) {
          cleaned[j] <- ""
          j <- j + 1
        }
        cleaned[j] <- line
        j <- j + 1
      } else {
        cleaned[j] <- remove_single_hash_comment(line)
        j <- j + 1
      }
    }
    cleaned <- cleaned[cleaned != ""]
    writeLines(cleaned, out_file)
  }
}


clean_project_code(input_folder, output_folder)


