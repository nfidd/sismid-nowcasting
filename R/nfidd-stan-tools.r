#' Find the Stan code in a package source tree
#'
#' Walks up from `path` looking for the source tree of this package, i.e. a
#' `DESCRIPTION` naming `nfidd.nowcasting` next to an `inst/stan` directory.
#' Returns `NULL` when there is no such directory above `path`.
#'
#' @param path Character string, the directory to start from.
#'
#' @return A character string with the path to `inst/stan`, or `NULL`.
#'
#' @keywords internal
.nfidd_source_stan_path <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = FALSE)
  repeat {
    description <- file.path(path, "DESCRIPTION")
    stan_path <- file.path(path, "inst", "stan")
    if (file.exists(description) && dir.exists(stan_path)) {
      package <- tryCatch(
        read.dcf(description, fields = "Package")[1, 1],
        error = function(e) NA_character_
      )
      if (isTRUE(package == "nfidd.nowcasting")) {
        return(stan_path)
      }
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      return(NULL)
    }
    path <- parent
  }
}

#' Get the path to Stan code
#'
#' This is the single place the Stan code is resolved, for both the models and
#' the functions they include. It looks in three places, in order.
#'
#' 1. The `nfidd.nowcasting.stan_path` option, if it is set. Set this to point
#'    the Stan tools at your own Stan code.
#' 2. The package source tree, if the working directory sits in a clone of the
#'    course repository. Rendering the course site, or running a session, from
#'    a clone therefore uses the Stan code in that clone rather than whichever
#'    version of the package happens to be installed. Without this a model
#'    added to `inst/stan` is invisible to [nfidd_stan_models()] and the rest
#'    of the Stan tools until the package is reinstalled. Using the source tree
#'    over an installed package is reported once per session, so it is never a
#'    silent switch. The report is only made in an interactive session, where
#'    someone is there to read it, and so never lands in a rendered page.
#' 3. The installed package.
#'
#' @return A character vector of directories holding the Stan code. This is a
#'   single directory unless the option names more than one.
#'
#' @family stantools
#'
#' @export
nfidd_stan_path <- function() {
  option_path <- getOption("nfidd.nowcasting.stan_path")
  if (!is.null(option_path)) {
    return(option_path)
  }

  source_path <- .nfidd_source_stan_path()
  installed_path <- system.file("stan", package = "nfidd.nowcasting")
  if (!is.null(source_path)) {
    if (nzchar(installed_path) && rlang::is_interactive()) {
      rlang::inform(
        c(
          paste0("Using the Stan code in ", source_path, "."),
          i = paste(
            "This is the course repository you are working in, not the",
            "installed package."
          ),
          i = paste(
            "Set the `nfidd.nowcasting.stan_path` option to use Stan code",
            "from somewhere else."
          )
        ),
        .frequency = "once",
        .frequency_id = "nfidd_stan_path_source_tree"
      )
    }
    return(source_path)
  }
  installed_path
}

#' Count the number of unmatched braces in a line
#' @noRd
.unmatched_braces <- function(line) {
  ifelse(
    grepl("{", line, fixed = TRUE),
    length(gregexpr("{", line, fixed = TRUE)), 0
  ) -
    ifelse(
      grepl("}", line, fixed = TRUE),
      length(gregexpr("}", line, fixed = TRUE)), 0
    )
}

#' Extract function names or content from Stan code
#'
#' @param content Character vector containing Stan code
#'
#' @param names_only Logical, if TRUE extract function names, otherwise
#'  extract function content.
#'
#' @param functions Optional, character vector of function names to extract
#'   content for.
#'
#' @return Character vector of function names or content
#'
#' @keywords internal
.extract_stan_functions <- function(
    content,
    names_only = FALSE,
    functions = NULL) {
  def_pattern <- "^(array\\[\\]\\s*)?(real|int|void|vector|row_vector|matrix)\\s+"
  func_pattern <- paste0(
    def_pattern,
    "(\\w+)\\s*\\("
  )
  func_lines <- grep(func_pattern, content, value = TRUE)
  # remove the func_pattern
  func_lines <- sub(def_pattern, "", func_lines)
  # get the next complete word after the pattern until the first (
  func_names <- sub("\\s*\\(.*$", "", func_lines)
  if (!is.null(functions)) {
    func_names <- intersect(func_names, functions)
  }
  if (names_only) {
    return(func_names)
  } else {
    func_content <- character(0)
    for (func_name in func_names) {
      start_line <- grep(paste0(def_pattern, func_name, "\\s*\\("), content)
      if (length(start_line) == 0) next
      end_line <- start_line
      brace_count <- 0
      # Ensure we find the first opening brace
      # Find first opening brace
      repeat {
        brace_count <- brace_count + .unmatched_braces(content[end_line])
        end_line <- end_line + 1
        if (brace_count > 0) break
      }

      # Continue until all braces are closed
      repeat {
        brace_count <- brace_count + .unmatched_braces(content[end_line])
        if (brace_count == 0) break
        end_line <- end_line + 1
      }

      func_content <- c(
        func_content,
        paste(content[start_line:end_line], collapse = "\n")
      )
    }
    return(func_content)
  }
}

#' Make Stan file paths relative to the Stan path they were found under
#'
#' `stan_path` can hold more than one directory, so each root is stripped in
#' turn. Matching is on a literal prefix rather than a pattern, so a root
#' holding regular expression characters cannot over-match.
#'
#' @param files Character vector of Stan file paths.
#' @param stan_path Character vector of Stan path roots.
#'
#' @return `files`, with any leading Stan path root removed.
#'
#' @keywords internal
.relative_to_stan_path <- function(files, stan_path) {
  for (root in stan_path) {
    prefix <- paste0(root, "/")
    under_root <- startsWith(files, prefix)
    files[under_root] <- substring(files[under_root], nchar(prefix) + 1)
  }
  files
}

#' Get Stan function names from Stan files
#'
#' This function reads all Stan files in the specified directory and extracts
#' the names of all functions defined in those files.
#'
#' @param stan_path Character vector of directories holding Stan files, as
#' returned by [nfidd_stan_path()]. May name more than one directory.
#'
#' @return A character vector containing unique names of all functions found in
#' the Stan files.
#'
#' @export
#'
#' @family stantools
nfidd_stan_functions <- function(
    stan_path = nfidd.nowcasting::nfidd_stan_path()) {
  stan_files <- list.files(
    file.path(stan_path, "functions"),
    pattern = "\\.stan$", full.names = TRUE,
    recursive = TRUE
  )
  functions <- character(0)
  for (file in stan_files) {
    content <- readLines(file)
    functions <- c(
      functions, .extract_stan_functions(content, names_only = TRUE)
    )
  }
  unique(functions)
}

#' Get Stan files containing specified functions
#'
#' This function retrieves Stan files from a specified directory, optionally
#' filtering for files that contain specific functions.
#'
#' @param functions Character vector of function names to search for. If NULL,
#'   all Stan files are returned.
#' @inheritParams nfidd_stan_functions
#'
#' @return A character vector of file paths to Stan files.
#'
#' @export
#'
#' @family stantools
nfidd_stan_function_files <- function(
    functions = NULL,
    stan_path = nfidd.nowcasting::nfidd_stan_path()) {
  # List all Stan files in the directory
  all_files <- list.files(
    file.path(stan_path, "functions"),
    pattern = "\\.stan$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (is.null(functions)) {
    return(all_files)
  } else {
    # Initialize an empty vector to store matching files
    matching_files <- character(0)

    for (file in all_files) {
      content <- readLines(file)
      extracted_functions <- .extract_stan_functions(content, names_only = TRUE)

      if (any(functions %in% extracted_functions)) {
        matching_files <- c(matching_files, file)
      }
    }

    # remove the Stan path from the file names
    matching_files <- .relative_to_stan_path(matching_files, stan_path)
    return(matching_files)
  }
}

#' Load Stan functions as a string
#'
#' @param functions Character vector of function names to load. Defaults to all
#' functions.
#'
#' @param stan_path Character vector of directories holding Stan files, as
#' returned by [nfidd_stan_path()]. May name more than one directory.
#'
#' @param wrap_in_block Logical, whether to wrap the functions in a
#' `functions{}` block. Default is FALSE.
#'
#' @param write_to_file Logical, whether to write the output to a file. Default
#' is FALSE.
#'
#' @param output_file Character string, the path to write the output file if
#' write_to_file is TRUE. Defaults to "nfidd_functions.stan".
#'
#' @return A character string containing the requested Stan functions
#'
#' @family stantools
#'
#' @export
nfidd_load_stan_functions <- function(
    functions = NULL, stan_path = nfidd.nowcasting::nfidd_stan_path(),
    wrap_in_block = FALSE, write_to_file = FALSE,
    output_file = "nfidd_functions.stan") {
  stan_files <- list.files(
    file.path(stan_path, "functions"),
    pattern = "\\.stan$", full.names = TRUE,
    recursive = TRUE
  )
  all_content <- character(0)

  for (file in stan_files) {
    content <- readLines(file)
    if (is.null(functions)) {
      all_content <- c(all_content, content)
    } else {
      for (func in functions) {
        func_content <- .extract_stan_functions(
          content,
          names_only = FALSE,
          functions = func
        )
        all_content <- c(all_content, func_content)
      }
    }
  }

  # Add version comment
  version_comment <- paste(
    "// Stan functions from nfidd version",
    utils::packageVersion("nfidd.nowcasting")
  )
  all_content <- c(version_comment, all_content)

  if (wrap_in_block) {
    all_content <- c("functions {", all_content, "}")
  }

  result <- paste(all_content, collapse = "\n")

  if (write_to_file) {
    writeLines(result, output_file)
    message("Stan functions written to: ", output_file, "\n")
  }

  return(result)
}

#' List Available Stan Models in NFIDD
#'
#' This function finds all available Stan models in the NFIDD package and
#' returns their names without the .stan extension.
#'
#' @param stan_path Character vector of directories holding Stan files, as
#'   returned by [nfidd_stan_path()]. May name more than one directory.
#'
#' @return A character vector of available Stan model names.
#'
#' @export
#'
#' @examples
#' nfidd_stan_models()
nfidd_stan_models <- function(stan_path = nfidd.nowcasting::nfidd_stan_path()) {
  stan_files <- list.files(
    stan_path,
    pattern = "\\.stan$", full.names = FALSE,
    recursive = FALSE
  )

  # Remove .stan extension. The same model can sit under more than one root.
  model_names <- unique(tools::file_path_sans_ext(stan_files))

  return(model_names)
}

#' Create a CmdStanModel with NFIDD Stan functions
#'
#' This function creates a CmdStanModel object using either a specified Stan
#' model from the NFIDD package or a custom Stan file provided by the user.
#'
#' @param model_name Character string specifying which Stan model to use from
#'   the NFIDD package. Ignored if model_file is provided.
#' @param model_file Character string specifying the path to a custom Stan file.
#'   If provided, this takes precedence over model_name.
#' @param include_paths Character vector of directories to include for Stan
#'   compilation. Defaults to the result of [nfidd_stan_path()], which also
#'   honours the "nfidd.nowcasting.stan_path" option.
#' @param ... Additional arguments passed to cmdstanr::cmdstan_model().
#'
#' @return A CmdStanModel object.
#'
#' @export
#'
#' @family modelhelpers
#'
#' @importFrom cmdstanr cmdstan_model
#'
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE)
#' if (!is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   # Using a model from the NFIDD package
#'   model <- nfidd_cmdstan_model("simple-nowcast", compile = FALSE)
#'   model
#'
#' }
nfidd_cmdstan_model <- function(
    model_name = NULL,
    model_file = NULL,
    include_paths = nfidd.nowcasting::nfidd_stan_path(),
    ...) {

  # Determine which Stan file to use
  if (!is.null(model_file)) {
    # Use custom model file
    if (!file.exists(model_file)) {
      stop(sprintf("Custom model file '%s' not found", model_file))
    }
    stan_model <- model_file
  } else if (!is.null(model_name)) {
    # Take the model from the same Stan code as the functions it includes. The
    # Stan path can hold more than one directory, so take the first one that
    # has the model in it.
    candidates <- file.path(
      nfidd.nowcasting::nfidd_stan_path(), paste0(model_name, ".stan")
    )
    found <- candidates[file.exists(candidates)]

    if (length(found) == 0) {
      stop(sprintf(
        "Model '%s.stan' not found in: %s",
        model_name,
        paste(nfidd.nowcasting::nfidd_stan_path(), collapse = ", ")
      ))
    }
    stan_model <- found[1]
  } else {
    stop("Either model_name or model_file must be provided")
  }

  cmdstan_model(
    stan_model,
    include_paths = include_paths,
    ...
  )
}

#' Sample from a CmdStanModel with NFIDD course defaults
#'
#' This function wraps the cmdstanr sample method with optimized defaults
#' for course use to speed up model fitting. All cmdstanr sample arguments
#' can still be overridden for experimentation.
#'
#' @param model A CmdStanModel object to sample from.
#' @param iter_warmup Integer, number of warmup iterations per chain.
#'   Defaults to 500 (reduced from cmdstanr default of 1000) for course speed.
#' @param iter_sampling Integer, number of sampling iterations per chain.
#'   Defaults to 500 (reduced from cmdstanr default of 1000) for course speed.
#' @param parallel_chains Integer, number of chains to run in parallel.
#'   Defaults to 4 for course speed.
#' @param save_warmup Logical, whether to save warmup samples.
#'   Defaults to FALSE for course speed.
#' @param ... Additional arguments passed to the model's sample method.
#'   All cmdstanr sample arguments are supported.
#'
#' @return A CmdStanMCMC object containing the posterior samples.
#'
#' @export
#'
#' @family modelhelpers
nfidd_sample <- function(model,
                         iter_warmup = 500,
                         iter_sampling = 500,
                         parallel_chains = 4,
                         save_warmup = FALSE,
                         ...) {
  model$sample(
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    parallel_chains = parallel_chains,
    save_warmup = save_warmup,
    ...
  )
}
