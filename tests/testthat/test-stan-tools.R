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
    expect_equal(nfidd_stan_models(), "a-model-that-is-not-installed")
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
