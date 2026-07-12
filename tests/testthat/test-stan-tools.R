#' Build a throwaway package source tree holding one Stan model
#' @noRd
local_source_tree <- function(model, package = "nfidd.nowcasting",
                              env = parent.frame()) {
  source_tree <- tempfile("source-tree")
  withr::defer(unlink(source_tree, recursive = TRUE), envir = env)
  stan_path <- file.path(source_tree, "inst", "stan")
  dir.create(stan_path, recursive = TRUE)
  writeLines(
    c(paste("Package:", package), "Version: 0.0.0.9000"),
    file.path(source_tree, "DESCRIPTION")
  )
  file.create(file.path(stan_path, paste0(model, ".stan")))
  source_tree
}

test_that(".nfidd_source_stan_path() finds the Stan code in a source tree", {
  source_tree <- local_source_tree("a-model")
  nested <- file.path(source_tree, "reference")
  dir.create(nested)

  expect_equal(
    normalizePath(.nfidd_source_stan_path(nested)),
    normalizePath(file.path(source_tree, "inst", "stan"))
  )
})

test_that("nfidd_stan_models() lists models the installed package lacks", {
  source_tree <- local_source_tree("a-model-that-is-not-installed")
  installed <- nfidd_stan_models(
    system.file("stan", package = "nfidd.nowcasting")
  )
  expect_false("a-model-that-is-not-installed" %in% installed)

  withr::with_dir(
    source_tree,
    expect_equal(
      suppressMessages(nfidd_stan_models()), "a-model-that-is-not-installed"
    )
  )
})

test_that("nfidd_stan_path() says when it uses a source tree", {
  source_tree <- local_source_tree("a-model")

  # defeat the once-per-session frequency so the test does not depend on
  # whether something else has already triggered the message
  withr::local_options(rlib_message_verbosity = "verbose")

  withr::with_options(
    list(rlang_interactive = TRUE),
    expect_message(
      withr::with_dir(source_tree, nfidd_stan_path()),
      "course repository"
    )
  )
})

test_that("nfidd_stan_path() stays quiet when nobody is there to read it", {
  source_tree <- local_source_tree("a-model")
  withr::local_options(rlib_message_verbosity = "verbose")

  # a rendered page is not an audience for the note
  withr::with_options(
    list(rlang_interactive = FALSE),
    expect_no_message(withr::with_dir(source_tree, nfidd_stan_path()))
  )
})

test_that("nfidd_stan_path() takes the option over everything else", {
  option_tree <- local_source_tree("an-option-model")
  source_tree <- local_source_tree("a-source-model")
  option_path <- file.path(option_tree, "inst", "stan")

  withr::with_options(
    list(nfidd.nowcasting.stan_path = option_path),
    withr::with_dir(source_tree, {
      expect_equal(nfidd_stan_path(), option_path)
      expect_equal(nfidd_stan_models(), "an-option-model")
    })
  )
})

test_that("nfidd_stan_path() falls back to the installed package", {
  outside <- tempfile("outside")
  on.exit(unlink(outside, recursive = TRUE), add = TRUE)
  dir.create(outside)

  withr::with_options(
    list(nfidd.nowcasting.stan_path = NULL),
    withr::with_dir(
      outside,
      expect_equal(
        nfidd_stan_path(),
        system.file("stan", package = "nfidd.nowcasting")
      )
    )
  )
})

test_that(".nfidd_source_stan_path() ignores other packages", {
  other <- local_source_tree("a-model", package = "someotherpackage")

  expect_null(.nfidd_source_stan_path(other))
})

test_that("nfidd_stan_models() lists the models used in the course", {
  models <- nfidd_stan_models()

  expect_true("estimate-inf-and-r-multi-stream" %in% models)
  expect_false(any(grepl("functions", models, fixed = TRUE)))
})

test_that("nfidd_stan_functions() includes functions used by the models", {
  functions <- nfidd_stan_functions()

  expect_true(all(
    c("renewal", "renewal_seeded", "geometric_random_walk") %in% functions
  ))
})

test_that("nfidd_stan_function_files() finds the file a function lives in", {
  expect_equal(
    nfidd_stan_function_files(functions = "renewal_seeded"),
    "functions/renewal_seeded.stan"
  )
})

test_that("nfidd_cmdstan_model() handles a Stan path with several entries", {
  extra <- local_source_tree("a-model")
  stan_path <- c(
    file.path(extra, "inst", "stan"),
    system.file("stan", package = "nfidd.nowcasting")
  )

  withr::with_options(list(nfidd.nowcasting.stan_path = stan_path), {
    expect_equal(nfidd_stan_path(), stan_path)

    model <- nfidd_cmdstan_model("simple-nowcast", compile = FALSE)
    expect_equal(
      normalizePath(model$stan_file()),
      normalizePath(file.path(stan_path[2], "simple-nowcast.stan"))
    )

    expect_error(
      nfidd_cmdstan_model("not-a-model", compile = FALSE),
      "not-a-model.stan' not found"
    )
  })
})

test_that("nfidd_stan_function_files() is relative under a multi-entry path", {
  extra <- local_source_tree("a-model")
  stan_path <- c(
    file.path(extra, "inst", "stan"),
    system.file("stan", package = "nfidd.nowcasting")
  )

  # renewal_seeded lives under the second directory only
  expect_equal(
    nfidd_stan_function_files(
      functions = "renewal_seeded", stan_path = stan_path
    ),
    "functions/renewal_seeded.stan"
  )
})

test_that(".relative_to_stan_path() strips roots literally", {
  # a root holding regular expression characters must not over-match
  expect_equal(
    .relative_to_stan_path(
      c("/a.b/stan/functions/renewal.stan", "/other/functions/renewal.stan"),
      c("/a.b/stan", "/other")
    ),
    c("functions/renewal.stan", "functions/renewal.stan")
  )
  expect_equal(
    .relative_to_stan_path("/axb/stan/f.stan", "/a.b/stan"),
    "/axb/stan/f.stan"
  )
})
