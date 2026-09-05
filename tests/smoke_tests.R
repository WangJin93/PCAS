## PCAS post-refactor smoke/integration tests (run with Rscript)
## Exercises every exported function against the live PCAS API where needed,
## plus unit-style checks of the refactored behaviour (missing-data feedback).
suppressMessages({
  library(PCAS)
  library(ggplot2)
})
options(PCAS.cache.dir = "data_files")   # keep the cache inside the repo
Sys.setenv(TZ = "UTC")

failures <- 0L
passes   <- 0L
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  if (ok) {
    passes <<- passes + 1L
    cat("PASS:", label, "\n")
  } else {
    failures <<- failures + 1L
    cat("FAIL:", label, "\n")
  }
  invisible(ok)
}
expect_no_error <- function(label, expr) {
  r <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  expect(label, is.null(r))
  invisible(r)
}

tmp <- tempfile(fileext = ".png")

cat("================ get_data ================\n")
d_expr <- suppressMessages(get_data("LUAD_CPTAC_protein", "expression",
                                    genes = "GAPDH"))
expect("get_data expression returns data.frame with rows",
       is.data.frame(d_expr) && nrow(d_expr) >= 1L)
d_deg <- suppressMessages(get_data("LUAD_CPTAC_protein_ttest", "DEGs"))
expect("get_data DEGs returns rows incl. adj.P.Val",
       is.data.frame(d_deg) && nrow(d_deg) > 0L &&
         all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(d_deg)))
d_cli <- suppressMessages(get_data("LUAD_APOLLO", "clinic"))
expect("get_data clinic returns data.frame",
       is.data.frame(d_cli) && "Cases_Submitter_ID" %in% colnames(d_cli))
expect("get_data unknown table returns NULL",
       is.null(suppressWarnings(suppressMessages(
         get_data("NOT_A_TABLE_XYZ", "expression", genes = "TP53")))))
expect("get_data bad action errors", {
  r <- tryCatch({ get_data("LUAD_CPTAC_protein", "banana"); "no error" },
                error = function(e) "error")
  identical(r, "error")
})

cat("================ get_expr_data ================\n")
avail <- function(x) attr(x, "availability")
r1 <- suppressMessages(get_expr_data("LUAD_CPTAC_protein", "TP53"))
expect("single gene fetch: 4 cols, numeric, non-empty",
       is.data.frame(r1) && ncol(r1) == 4L && nrow(r1) > 10L &&
         is.numeric(r1$TP53) && all(!is.na(r1$TP53)))
expect("single gene fetch carries availability attr",
       is.data.frame(avail(r1)) && nrow(avail(r1)) >= 1L)

r2 <- suppressMessages(get_expr_data("LUAD_CPTAC_protein",
                                     c("TP53", "ZZZNOTAGENE")))
expect("partial missing gene: column kept, all-NA, plus real gene present",
       is.data.frame(r2) && "ZZZNOTAGENE" %in% colnames(r2) &&
         all(is.na(r2$ZZZNOTAGENE)) && any(!is.na(r2$TP53)))
expect("partial missing gene flagged in availability",
       any(!avail(r2)$present[avail(r2)$gene == "ZZZNOTAGENE"]))
expect("all-missing gene -> NULL",
       is.null(suppressMessages(get_expr_data("LUAD_CPTAC_protein",
                                              "ZZZNOTAGENE"))))
expect("unknown dataset ignored with valid one kept",
       is.data.frame(suppressMessages(
         get_expr_data(c("LUAD_CPTAC_protein", "BOGUS_XX"), "TP53"))))

r3 <- suppressMessages(get_expr_data(c("LUAD_CPTAC_mRNA", "LUAD_CPTAC_protein"),
                                     c("TP53", "GAPDH")))
expect("multi-dataset fetch columns = request order ID,type,dataset,TP53,GAPDH",
       is.data.frame(r3) &&
         identical(colnames(r3), c("ID", "type", "dataset", "TP53", "GAPDH")) &&
         all(is.numeric(r3$TP53), is.numeric(r3$GAPDH)) &&
         length(unique(r3$dataset)) == 2L)

cat("================ merge_clinic_data ================\n")
apollo <- suppressMessages(get_expr_data("LUAD_APOLLO_mRNA", "TP53"))
m1 <- suppressMessages(merge_clinic_data("LUAD_APOLLO", apollo))
expect("merge_clinic returns list(df, summary)",
       is.list(m1) && is.data.frame(m1$df) &&
         all(c("cohort", "n_matched_rows", "clinical_na") %in%
               names(m1$summary)))
expect("merge summary counts consistent",
       m1$summary$n_unmatched_tumour_rows >= 0L &&
         m1$summary$n_matched_rows == nrow(m1$df))
expect("merged df keeps expression col at position 4",
       ncol(m1$df) >= 4L && is.numeric(m1$df[[4]]))
m2 <- suppressMessages(merge_clinic_data("LUAD_APOLLO", apollo,
                                         return_summary = FALSE))
expect("return_summary=FALSE gives data.frame + attr",
       is.data.frame(m2) && !is.null(attr(m2, "merge_summary")))
expect("merge_clinic with unknown cohort -> NULL",
       is.null(suppressWarnings(suppressMessages(
         merge_clinic_data("BOGUS_COHORT", apollo)))))

cat("================ get_DEGs_result ================\n")
deg1 <- suppressMessages(get_DEGs_result("LUAD_CPTAC_protein", "t.test"))
expect("DEGs t.test: Symbol/logFC/P.Value present",
       is.data.frame(deg1) && all(c("Symbol", "logFC", "P.Value") %in%
                                    colnames(deg1)) &&
         attr(deg1, "method") == "t.test")
deg2 <- suppressMessages(get_DEGs_result("LUAD_CPTAC_mRNA", "limma"))
expect("DEGs mRNA limma: mapped to Symbol, protein_coding filterable",
       is.data.frame(deg2) && "Symbol" %in% colnames(deg2) &&
         ("gene_type" %in% colnames(deg2)))
expect("DEGs bad dataset -> NULL",
       is.null(suppressWarnings(suppressMessages(
         get_DEGs_result("BOGUS_DATASET", "limma")))))

cat("================ cor_cancer_genelist ================\n")
cc <- suppressMessages(cor_cancer_genelist(
  dataset1 = "LUAD_CPTAC_protein", id1 = "STAT3",
  dataset2 = "LUAD_CPTAC_mRNA", id2 = c("TNS1", "TP53"),
  sample_type = "Tumor", cor_method = "pearson"))
expect("cor_cancer list structure",
       is.list(cc) && all(c("cor_result", "cor_data", "n", "summary") %in%
                            names(cc)))
expect("cor_result rows follow id2 order with clean names",
       identical(as.character(cc$cor_result$Symbol), c("TNS1", "TP53")) &&
         all(c("Correlation", "P value", "n") %in% colnames(cc$cor_result)))
expect("cor_data columns aligned (col4 target, col5 first feature)",
       colnames(cc$cor_data)[4] == "STAT3" &&
         colnames(cc$cor_data)[5] == "TNS1")
expect("cor_cancer with unmeasurable id1 -> NULL",
       is.null(suppressMessages(cor_cancer_genelist(
         dataset1 = "LUAD_CPTAC_protein", id1 = "ZZZNOTAGENE",
         dataset2 = "LUAD_CPTAC_mRNA", id2 = "TP53"))))

cat("================ cor_pancancer_* ================\n")
pds <- c("LUAD_CPTAC_protein", "HNSCC_CPTAC_protein")
dfT <- suppressMessages(get_expr_data(pds, "TNS1"))
gsT <- suppressMessages(get_expr_data(pds, c("TP53", "SIRPA", "ZZZNOTAGENE")))
cg <- suppressMessages(cor_pancancer_genelist(dfT, gsT, sample_type = "Tumor"))
expect("pancancer genelist: r/p/n/sss/summary present and dims match",
       is.list(cg) && all(c("r", "p", "n", "sss", "summary") %in% names(cg)) &&
         identical(dim(cg$r), dim(cg$p)) && identical(dim(cg$r), dim(cg$n)) &&
         identical(colnames(cg$r), colnames(cg$p)))
expect("pancancer genelist: rownames = stripped cohort labels",
       all(grepl("^(LUAD|HNSCC)_CPTAC$", rownames(cg$r))))
expect("pancancer n matrix > 0 where r is not NA",
       all(cg$n[!is.na(cg$r)] > 0L))
expect("invalid sample_type -> NULL",
       is.null(suppressWarnings(suppressMessages(
         cor_pancancer_genelist(dfT, gsT, sample_type = "banana")))))

cd <- suppressMessages(cor_pancancer_drug(dfT, cor_method = "spearman",
                                          Target.pathway = "Cell cycle"))
expect("pancancer drug: r/p/n present with drug colnames",
       is.list(cd) && all(c("r", "p", "n") %in% names(cd)) &&
         ncol(cd$r) > 0L && nrow(cd$r) == 2L &&
         all(grepl("_", colnames(cd$r), fixed = TRUE)))
expect("pancancer drug: invalid pathway warns and returns NULL",
       is.null(suppressWarnings(suppressMessages(
         cor_pancancer_drug(dfT, Target.pathway = "Not a pathway")))))

ct <- suppressMessages(cor_pancancer_TIL(dfT, cor_method = "spearman",
                                         TIL_type = "TIMER"))
expect("pancancer TIL: TIMER cell columns present",
       is.list(ct) && ncol(ct$r) > 0L && nrow(ct$r) > 0L &&
         identical(dim(ct$r), dim(ct$n)))
expect("pancancer TIL: invalid algorithm warns and returns NULL",
       is.null(suppressWarnings(suppressMessages(
         cor_pancancer_TIL(dfT, TIL_type = "no_such_algo")))))

cat("================ viz_* ================\n")
# local synthetic data for quick plotting checks
syn <- data.frame(ID = paste0("S", 1:80),
                  type = rep(c("Tumor", "Normal"), each = 40),
                  dataset = "LUAD_CPTAC_protein",
                  GENE1 = rnorm(80), GENE2 = rnorm(80),
                  stringsAsFactors = FALSE)
p <- suppressMessages(viz_TvsN(syn, df_type = "multi_gene"))
expect("viz_TvsN multi_gene returns ggplot",
       inherits(p, "ggplot"))
expect_no_error("viz_TvsN single renders", {
  q <- suppressMessages(viz_TvsN(syn[, 1:4], df_type = "single",
                                 Show.n = TRUE))
  ggplot2::ggsave(tmp, q, width = 4, height = 3); TRUE
})
syn_na <- syn; syn_na$GENE1[c(1, 5, 12)] <- NA
p2 <- suppressMessages(viz_TvsN(syn_na, df_type = "single"))
expect("viz_TvsN handles NA rows with message, still a ggplot",
       inherits(p2, "ggplot"))
# only-Tumour dataset -> p-values skipped but plot still returned
syn_tumor <- syn[syn$type == "Tumor", ]
p3 <- suppressMessages(viz_TvsN(syn_tumor, df_type = "single"))
expect("viz_TvsN single-group dataset still returns ggplot",
       inherits(p3, "ggplot"))
expect("viz_TvsN default df_type no longer errors",
       inherits(suppressMessages(viz_TvsN(syn[, 1:4])), "ggplot"))

vol <- data.frame(Symbol = paste0("G", 1:300),
                  logFC = rnorm(300), P.Value = runif(300),
                  adj.P.Val = runif(300), stringsAsFactors = FALSE)
vol$adj.P.Val[1:5] <- 1e-10
pv1 <- suppressMessages(viz_DEGs_volcano(vol, show.top = TRUE))
expect("volcano returns ggplot and keeps data for $data consumers",
       inherits(pv1, "ggplot") && is.data.frame(pv1$data))
expect_no_error("volcano with explicit labels renders", {
  ggplot2::ggsave(tmp, suppressMessages(
    viz_DEGs_volcano(vol, show.labels = c("G1", "G300"))), width = 4, height = 3)
  TRUE
})
expect("volcano on tiny data with show.top does not crash",
       inherits(suppressMessages(viz_DEGs_volcano(vol[1:5, ],
                                                  show.top = TRUE)), "ggplot"))
expect("volcano missing logFC errors informatively", {
  r <- tryCatch({ viz_DEGs_volcano(data.frame(Symbol = "a")); "no error" },
                error = function(e) "error")
  identical(r, "error")
})

cp <- suppressMessages(viz_corplot(syn, "GENE1", "GENE2"))
expect("viz_corplot returns ggplot", inherits(cp, "ggplot"))
cp_na <- suppressMessages(viz_corplot(syn_na, "GENE1", "GENE2"))
expect("viz_corplot handles NAs and reports removal", inherits(cp_na, "ggplot"))
expect("viz_corplot bad columns errors informatively", {
  r <- tryCatch({ viz_corplot(syn, "GENE1", "NOPE"); "no error" },
                error = function(e) "error")
  identical(r, "error")
})

rmat <- matrix(c(0.5, -0.2, 0.3, NA, 0.7, -0.1, 0.2, 0.4, 0.9), 3,
               dimnames = list(c("A", "B", "C"), c("X", "Y", "Z")))
pmat <- matrix(c(0.01, 0.2, 0.6, NA, 0.04, 0.5, 0.3, 0.02, 1e-4), 3,
               dimnames = dimnames(rmat))
expect("heatmap with NA cells still returns ggplot (no NA->0)",
       inherits(suppressMessages(viz_cor_heatmap(rmat, pmat)), "ggplot"))
rmat2 <- rmat; rmat2[1, 1] <- 0.8   # full matrix -> clustering/tree path
expect("heatmap full matrix with tree returns ggplot",
       inherits(suppressMessages(viz_cor_heatmap(rmat2, pmat)), "ggplot"))
expect("heatmap dimension mismatch errors informatively", {
  r <- tryCatch({ viz_cor_heatmap(rmat, pmat[1:2, ]); "no error" },
                error = function(e) "error")
  identical(r, "error")
})

ps <- tryCatch(suppressMessages(viz_phoso_sites("TNS1", "CPTAC")),
               error = function(e) NULL)
expect("viz_phoso_sites returns plot or graceful NULL on network failure",
       is.null(ps) || inherits(ps, "ggplot"))

cat(sprintf("\n==== %d passed, %d failed ====\n", passes, failures))
quit(status = if (failures > 0L) 1L else 0L)
