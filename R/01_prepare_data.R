suppressPackageStartupMessages(library(data.table))
source("R/utils.R")
ensure_directories()

label_path <- "data/raw/deepsynergy_labels.csv"
drugcomb_path <- "data/raw/drugcomb_v15_oneil.csv"
identifier_path <- "data/interim/drugcomb_oneil_identifiers.csv"

for (path in c(label_path, drugcomb_path, identifier_path)) {
  assert_true(file.exists(path), paste("Required source file is missing:", path))
}

labels <- fread(label_path)
if (names(labels)[1L] == "") {
  setnames(labels, 1L, "source_triplet_label")
}
required_label_columns <- c("drug_a_name", "drug_b_name", "cell_line", "synergy", "fold")
assert_true(all(required_label_columns %in% names(labels)), "DeepSynergy labels are missing required columns")
assert_true(nrow(labels) == 23052L, "Unexpected DeepSynergy label row count")
assert_true(all(is.finite(labels$synergy)), "DeepSynergy labels contain non-finite synergy values")

labels[, pair_id := canonical_pair_id(drug_a_name, drug_b_name)]
labels[, triplet_id := triplet_id(drug_a_name, drug_b_name, cell_line)]

pair_fold_check <- labels[, .(fold_count = uniqueN(fold)), by = pair_id]
assert_true(all(pair_fold_check$fold_count == 1L), "At least one drug pair occurs in multiple outer folds")
assert_true(setequal(unique(labels$fold), 0:4), "Published folds are not numbered 0 through 4")
assert_true(uniqueN(labels$pair_id) == 583L, "Unexpected number of drug pairs")
assert_true(uniqueN(c(labels$drug_a_name, labels$drug_b_name)) == 38L, "Unexpected number of drugs")
assert_true(uniqueN(labels$cell_line) == 39L, "Unexpected number of cell lines")
assert_true(uniqueN(labels$triplet_id) == 22737L, "Unexpected number of unique triplets")

outcomes <- labels[, .(
  drug_1 = min(c(drug_a_name, drug_b_name)),
  drug_2 = max(c(drug_a_name, drug_b_name)),
  cell_line = first(cell_line),
  loewe = mean(synergy),
  loewe_sd = finite_sd(synergy),
  loewe_min = min(synergy),
  loewe_max = max(synergy),
  loewe_range = max(synergy) - min(synergy),
  replicate_count = .N,
  outer_fold = first(fold)
), by = .(triplet_id, pair_id)]

assert_true(nrow(outcomes) == 22737L, "Triplet aggregation did not produce 22,737 rows")
assert_true(sum(outcomes$replicate_count == 2L) == 315L, "Unexpected number of repeated triplets")
assert_true(all(outcomes$replicate_count %in% c(1L, 2L)), "Unexpected replicate count")

drug_map <- data.table(
  deepsynergy_name = c(
    "5-FU", "ABT-888", "AZD1775", "BEZ-235", "BORTEZOMIB", "CARBOPLATIN",
    "CYCLOPHOSPHAMIDE", "DASATINIB", "DEXAMETHASONE", "DINACICLIB", "DOXORUBICIN",
    "ERLOTINIB", "ETOPOSIDE", "GELDANAMYCIN", "GEMCITABINE", "L778123", "LAPATINIB",
    "METFORMIN", "METHOTREXATE", "MITOMYCINE", "MK-2206", "MK-4541", "MK-4827",
    "MK-5108", "MK-8669", "MK-8776", "MRK-003", "OXALIPLATIN", "PACLITAXEL",
    "PD325901", "SN-38", "SORAFENIB", "SUNITINIB", "TEMOZOLOMIDE", "TOPOTECAN",
    "VINBLASTINE", "VINORELBINE", "ZOLINZA"
  ),
  drugcomb_name = c(
    "5-Fluorouracil", "Veliparib", "MK-1775", "915019-65-7", "Bortezomib", "Carboplatinum",
    "cyclophosphamide", "Dasatinib", "dexamethasone", "Dinaciclib", "doxorubicin",
    "Erlotinib", "etoposide", "geldanamycin", "gemcitabine", "L-778123 free base", "Lapatinib",
    "metformin", "methotrexate", "mitomycin C", "MK-2206", "MK-4541", "Niraparib",
    "MK-5108", "Deforolimus", "891494-63-6", "MRK003", "Eloxatin (TN)", "paclitaxel",
    "391210-10-9", "7-Ethyl-10-hydroxycamptothecin", "Sorafenib", "Sunitinib", "temozolomide",
    "topotecan", "vinblastine", "vinorelbine", "Vorinostat"
  )
)
drug_map[, mapping_basis := "exact or official DrugComb synonym"]
drug_map[deepsynergy_name %in% c("CARBOPLATIN", "MK-8776", "OXALIPLATIN"),
  mapping_basis := "manual check against official DrugComb identifier row"]

identifier_table <- fread(identifier_path)
assert_true(uniqueN(identifier_table$dname) == 38L, "Identifier table does not contain 38 unique drugs")
drug_map <- merge(
  drug_map,
  identifier_table[, .(
    drugcomb_name = dname,
    drugcomb_id = id,
    pubchem_cid = cid,
    inchikey,
    chembl_id,
    drugbank_id,
    clinical_phase
  )],
  by = "drugcomb_name",
  all.x = TRUE,
  sort = FALSE
)
assert_true(!anyNA(drug_map$drugcomb_id), "At least one drug map lacks an official DrugComb identifier")
assert_true(setequal(drug_map$deepsynergy_name, unique(c(labels$drug_a_name, labels$drug_b_name))),
  "Drug map does not match all DeepSynergy names")

drugcomb <- fread(drugcomb_path, na.strings = c("NA", "\\N"))
required_drugcomb_columns <- c(
  "block_id", "drug_row", "drug_col", "cell_line_name", "study_name",
  "css_ri", "synergy_loewe"
)
assert_true(all(required_drugcomb_columns %in% names(drugcomb)), "DrugComb subset is missing required columns")
assert_true(nrow(drugcomb) == 92208L, "Unexpected DrugComb O'Neil row count")
assert_true(all(drugcomb$study_name == "ONEIL"), "DrugComb subset contains another study")

reverse_drug_map <- setNames(drug_map$deepsynergy_name, drug_map$drugcomb_name)
drugcomb[, drug_1_deepsynergy := unname(reverse_drug_map[drug_row])]
drugcomb[, drug_2_deepsynergy := unname(reverse_drug_map[drug_col])]
assert_true(!anyNA(drugcomb$drug_1_deepsynergy) && !anyNA(drugcomb$drug_2_deepsynergy),
  "At least one DrugComb drug name was not mapped")

cell_names <- sort(unique(labels$cell_line))
cell_keys <- normalize_identifier(cell_names)
assert_true(!anyDuplicated(cell_keys), "Normalized DeepSynergy cell-line names are not unique")
cell_map <- setNames(cell_names, cell_keys)
drugcomb[, cell_line := unname(cell_map[normalize_identifier(cell_line_name)])]
assert_true(!anyNA(drugcomb$cell_line), "At least one DrugComb cell-line name was not mapped")

drugcomb[, pair_id := canonical_pair_id(drug_1_deepsynergy, drug_2_deepsynergy)]
drugcomb[, triplet_id := triplet_id(drug_1_deepsynergy, drug_2_deepsynergy, cell_line)]

drugcomb_triplets <- drugcomb[, .(
  drugcomb_block_count = .N,
  css = mean(css_ri),
  css_sd = finite_sd(css_ri),
  css_min = min(css_ri),
  css_max = max(css_ri),
  css_range = max(css_ri) - min(css_ri),
  drugcomb_loewe = mean(synergy_loewe),
  drugcomb_loewe_sd = finite_sd(synergy_loewe)
), by = .(triplet_id, pair_id)]

assert_true(nrow(drugcomb_triplets) == 22737L, "DrugComb aggregation did not produce 22,737 triplets")
assert_true(all(drugcomb_triplets$drugcomb_block_count %in% c(4L, 8L)),
  "DrugComb triplets do not have four blocks per source replicate")

analysis_data <- merge(
  outcomes,
  drugcomb_triplets,
  by = c("triplet_id", "pair_id"),
  all = TRUE,
  sort = FALSE
)
assert_true(nrow(analysis_data) == 22737L, "DeepSynergy and DrugComb triplets did not join one to one")
assert_true(!anyNA(analysis_data$loewe) && !anyNA(analysis_data$css),
  "At least one outcome was unmatched")
assert_true(all(analysis_data$drugcomb_block_count == 4L * analysis_data$replicate_count),
  "DrugComb block counts do not agree with DeepSynergy replicate counts")

setorder(analysis_data, outer_fold, pair_id, cell_line)
fwrite(analysis_data, "data/processed/oneil_analysis_outcomes.csv")
saveRDS(analysis_data, "data/processed/oneil_analysis_outcomes.rds", compress = "xz")
fwrite(drug_map[order(deepsynergy_name)], "outputs/tables/drug_name_mapping.csv")

cell_mapping <- unique(drugcomb[, .(
  deepsynergy_cell_line = cell_line,
  drugcomb_cell_line = cell_line_name
)])
setorder(cell_mapping, deepsynergy_cell_line)
fwrite(cell_mapping, "outputs/tables/cell_line_mapping.csv")

duplicate_audit <- analysis_data[replicate_count > 1L, .(
  drug_1, drug_2, cell_line, pair_id, outer_fold,
  replicate_count, loewe, loewe_sd, loewe_min, loewe_max, loewe_range,
  drugcomb_block_count, css, css_sd, css_range
)]
fwrite(duplicate_audit, "outputs/tables/duplicate_triplet_audit.csv")

audit_summary <- data.table(
  item = c(
    "DeepSynergy label rows",
    "Unique drug-pair-cell-line triplets",
    "Triplets with two source measurements",
    "Unique drug pairs",
    "Unique drugs",
    "Unique cell lines",
    "Pairs assigned to more than one outer fold",
    "DrugComb v1.5 O'Neil block rows",
    "DrugComb triplets matched to DeepSynergy",
    "Triplets with missing CSS",
    "Pearson correlation: published and DrugComb Loewe",
    "Spearman correlation: published and DrugComb Loewe",
    "Median range between repeated published Loewe scores",
    "Maximum range between repeated published Loewe scores"
  ),
  value = c(
    nrow(labels),
    nrow(analysis_data),
    sum(analysis_data$replicate_count == 2L),
    uniqueN(analysis_data$pair_id),
    uniqueN(c(analysis_data$drug_1, analysis_data$drug_2)),
    uniqueN(analysis_data$cell_line),
    sum(pair_fold_check$fold_count > 1L),
    nrow(drugcomb),
    sum(!is.na(analysis_data$css)),
    sum(is.na(analysis_data$css)),
    cor(analysis_data$loewe, analysis_data$drugcomb_loewe, method = "pearson"),
    cor(analysis_data$loewe, analysis_data$drugcomb_loewe, method = "spearman"),
    median(analysis_data[replicate_count == 2L]$loewe_range),
    max(analysis_data[replicate_count == 2L]$loewe_range)
  )
)
fwrite(audit_summary, "outputs/tables/data_audit_summary.csv")

source_manifest <- data.table(
  source = c(
    "DeepSynergy labels",
    "DeepSynergy published feature matrix",
    "DrugComb v1.5 full summary archive",
    "DrugComb v1.5 filtered O'Neil subset",
    "DrugComb v1.4 drug identifiers"
  ),
  revision = c(
    "tastanlab/ohe commit 4d3368d0c004949e319842c56f7900fe3ccd5f4b",
    "Hugging Face commit 2cc715df1e91d0a42269cb2f58ab32fce30ad2dd",
    "Zenodo record 15235991, version 1.5",
    "Derived by exact study_name filter from Zenodo record 15235991",
    "Zenodo record 18449193, version 1.4"
  ),
  local_file = c(
    label_path,
    "data/raw/deepsynergy_X.p.gz",
    NA_character_,
    drugcomb_path,
    "data/raw/DrugComb_drug_identifiers_v1.4.xlsx"
  ),
  checksum = c(
    "sha256:4ef40678fff51faa84fa4de63344c5b51e2e4abe7c28a417af24ec2b50fa135a",
    "sha256:e399f6cdb40184048ac34b370f89cc6e4582bec72893cc2c98d90c1634517899",
    "md5:cd3e24758c10ffaea0124b6e80626c32",
    "sha256:6bf6e7e434e4866d48c4e455f5a1039ed4bc4cf355c1f8022bcef71c9ea96349",
    "md5:06894e63e5e06baa23d2bd59d6bd8808"
  ),
  url = c(
    "https://raw.githubusercontent.com/tastanlab/ohe/4d3368d0c004949e319842c56f7900fe3ccd5f4b/DeepSynergy/labels.csv",
    "https://huggingface.co/datasets/ebcandir/DeepSynergy/resolve/2cc715df1e91d0a42269cb2f58ab32fce30ad2dd/X.p.gz",
    "https://zenodo.org/records/15235991",
    "https://zenodo.org/api/records/15235991/files/summary_v_1_5.csv/content",
    "https://zenodo.org/api/records/18449193/files/DrugComb_drug_identifiers.xlsx/content"
  )
)
fwrite(source_manifest, "data/raw/SOURCE_MANIFEST.csv")

log_lines <- c(
  sprintf("Prepared %s unique triplets from %s label rows.", nrow(analysis_data), nrow(labels)),
  sprintf("Aggregated %s repeated triplets.", sum(analysis_data$replicate_count == 2L)),
  sprintf("Verified %s pair-disjoint outer folds across %s pairs.", uniqueN(analysis_data$outer_fold), uniqueN(analysis_data$pair_id)),
  sprintf("Matched CSS for %s triplets from %s DrugComb block rows.", sum(!is.na(analysis_data$css)), nrow(drugcomb)),
  sprintf("Published versus DrugComb Loewe Pearson correlation: %.4f.", cor(analysis_data$loewe, analysis_data$drugcomb_loewe))
)
writeLines(log_lines, "logs/01_prepare_data.log")
cat(paste(log_lines, collapse = "\n"), "\n")
