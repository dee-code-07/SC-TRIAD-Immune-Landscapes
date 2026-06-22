# sc-triad patch: harmonize cell_type labels in triad_integrated.rds
# deeksha h | reg. 241706005 | msc bioinformatics iv sem
#
# problem: script 04 created cell_type_harmonized internally but never wrote it
#          back to the rds. cell_type in the saved object still carries:
#            "Platelet"  -- only in t2d/control origin cells
#            "CD16 NK"   -- only in t2d/control origin cells
#          this caused cellchat (script 06) to treat Platelet and Megakaryocyte
#          as distinct cell types across conditions, creating a false signal of
#          complete cell type loss in t2d and gain in asthma.
#
# fix: collapse Platelet -> Megakaryocyte, CD16 NK -> NK in cell_type column.
#      overwrite triad_integrated.rds in-place.
#      scripts 01-05 are unaffected. script 06 must be re-run after this patch.
#
# runtime: < 2 minutes (no computation, just label assignment + saveRDS)

suppressPackageStartupMessages(library(Seurat))

rds_path <- "/home/deekshah/sc-triad/02_scrna/pbmc/03_integration/triad_integrated.rds"

cat(format(Sys.time(), "[%H:%M:%S]"), "loading integrated object...\n")
obj <- readRDS(rds_path)
cat(format(Sys.time(), "[%H:%M:%S]"), "cells:", ncol(obj), "\n")

cat("\ncell_type distribution before patch:\n")
print(sort(table(obj$cell_type), decreasing = TRUE))

# verify the labels that need patching are present
before_platelet <- sum(obj$cell_type == "Platelet",  na.rm = TRUE)
before_cd16nk   <- sum(obj$cell_type == "CD16 NK",   na.rm = TRUE)
cat("\nPlatelet cells to relabel:", before_platelet, "\n")
cat("CD16 NK cells to relabel:", before_cd16nk, "\n")

# apply harmonization
obj$cell_type[obj$cell_type == "Platelet"]  <- "Megakaryocyte"
obj$cell_type[obj$cell_type == "CD16 NK"]   <- "NK"

cat("\ncell_type distribution after patch:\n")
print(sort(table(obj$cell_type), decreasing = TRUE))

# verify zero residual
stopifnot(sum(obj$cell_type == "Platelet", na.rm = TRUE) == 0)
stopifnot(sum(obj$cell_type == "CD16 NK",  na.rm = TRUE) == 0)
cat("\nverification passed: no Platelet or CD16 NK labels remain\n")

cat(format(Sys.time(), "[%H:%M:%S]"), "saving patched rds...\n")
saveRDS(obj, rds_path)
cat(format(Sys.time(), "[%H:%M:%S]"), "done:", rds_path, "\n")
cat("next step: run 06_cellchat.R\n")
