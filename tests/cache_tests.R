## PCAS cache behaviour tests (GCAS-style local caching)
## Verifies that an identical repeated request is served from the local cache
## and performs ZERO network calls.
suppressMessages(library(PCAS))

failures <- 0L; passes <- 0L
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  if (ok) { passes <<- passes + 1L; cat("PASS:", label, "\n") }
  else    { failures <<- failures + 1L; cat("FAIL:", label, "\n") }
  invisible(ok)
}

## ---- helper: block the network layer inside the PCAS namespace ----------
block_get_data <- function(counter) {
  ns <- asNamespace("PCAS")
  old <- get("get_data", envir = ns)
  unlockBinding("get_data", ns)
  mock <- function(...) {
    counter$n <- counter$n + 1L
    stop("get_data() was called although the data should come from the cache")
  }
  assign("get_data", mock, envir = ns)
  lockBinding("get_data", ns)
  invisible(old)
}
restore_get_data <- function(old) {
  ns <- asNamespace("PCAS")
  unlockBinding("get_data", ns)
  assign("get_data", old, envir = ns)
  lockBinding("get_data", ns)
}

tmp <- tempfile("pcas_cache_test_")
dir.create(tmp, recursive = TRUE)

## ================= 1. expression cache (get_expr_data) ===================
cat("--------------- get_expr_data cache ---------------\n")
r1 <- suppressMessages(get_expr_data("LUAD_CPTAC_protein", c("TP53", "TNS1"),
                                     cache_dir = tmp))
expect("first expression query returns data",
       is.data.frame(r1) && nrow(r1) > 10L)
expect("expression cache file was created",
       length(list.files(file.path(tmp, "data_temp"),
                         pattern = "\\.RData$")) >= 1L)

counter <- new.env(); counter$n <- 0L
old <- block_get_data(counter)
r2 <- suppressMessages(get_expr_data("LUAD_CPTAC_protein", c("TP53", "TNS1"),
                                     cache_dir = tmp))
restore_get_data(old)
expect("identical repeated expression query is served from cache (no network)",
       is.data.frame(r2) && nrow(r2) == nrow(r1) && counter$n == 0L)
expect("cached and freshly-fetched expression results are equal",
       isTRUE(all.equal(r1, r2)))

## an uncached gene combination must still hit the network (counter + 1) and
## fail when the network is blocked -> proves it really tried to fetch
counter <- new.env(); counter$n <- 0L
old <- block_get_data(counter)
r3 <- tryCatch(suppressMessages(
  get_expr_data("LUAD_CPTAC_protein", "TP53", cache_dir = tmp)),
  error = function(e) NULL)
restore_get_data(old)
expect("different (uncached) query attempts the network exactly once",
       is.null(r3) && counter$n == 1L)

## ================= 2. DEGs cache (get_DEGs_result) ========================
cat("--------------- get_DEGs_result cache ---------------\n")
d1 <- suppressMessages(get_DEGs_result("LUAD_CPTAC_protein", "t.test",
                                       cache_dir = tmp))
expect("first DEGs query returns data",
       is.data.frame(d1) && nrow(d1) > 100L)
expect("DEGs cache file was created",
       length(list.files(file.path(tmp, "DEG_results"),
                         pattern = "\\.RData$")) >= 1L)
counter <- new.env(); counter$n <- 0L
old <- block_get_data(counter)
d2 <- suppressMessages(get_DEGs_result("LUAD_CPTAC_protein", "t.test",
                                       cache_dir = tmp))
restore_get_data(old)
expect("identical repeated DEGs query is served from cache (no network)",
       is.data.frame(d2) && nrow(d2) == nrow(d1) && counter$n == 0L)
expect("cached and freshly-fetched DEGs results are equal",
       isTRUE(all.equal(d1, d2)))

## ================= 3. clinical table cache (merge_clinic_data) ============
cat("--------------- merge_clinic_data cache ---------------\n")
expr_fake <- data.frame(
  ID = paste0("C3L-", sprintf("%05d", 1:30), "_FAKE"),
  type = "Tumor", dataset = "LUAD_APOLLO",
  TP53 = rnorm(30), stringsAsFactors = FALSE)
m1 <- suppressMessages(suppressWarnings(
  merge_clinic_data("LUAD_APOLLO", expr_fake, cache_dir = tmp)))
expect("first clinical merge returns list(df, summary)",
       is.list(m1) && all(c("df", "summary") %in% names(m1)))
expect("clinical cache file was created",
       file.exists(file.path(tmp, "clinic_data", "LUAD_APOLLO.rds")))

counter <- new.env(); counter$n <- 0L
old <- block_get_data(counter)
m2 <- suppressMessages(suppressWarnings(
  merge_clinic_data("LUAD_APOLLO", expr_fake, cache_dir = tmp)))
restore_get_data(old)
expect("repeated clinical merge is served from cache (no network)",
       is.list(m2) && counter$n == 0L)

cat(sprintf("\n==== cache tests: %d passed, %d failed ====\n",
            passes, failures))
quit(status = if (failures > 0L) 1L else 0L)
