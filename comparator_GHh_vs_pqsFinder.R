setwd('/home/dimitri/rishav/table_comparator_G4hunt_pqsfinder')
getwd()

library(devtools)

# --------------------------------------------------------------------
# # Setup & required packages
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")

# # # Install if missing
# pkgs <- c("readxl","dplyr","GenomicRanges","BSgenome.Hsapiens.UCSC.hg38",
#           "Biostrings","pqsfinder")
# for (p in pkgs) if (!requireNamespace(p, quietly=TRUE))
#   BiocManager::install(p)
# 
# devtools::install_github("rongxinzh/G4SNVHunter", dependencies = TRUE)

# 
# if (!requireNamespace("pwalign", quietly=TRUE)) {
#   BiocManager::install("pwalign")
# }

library(readxl)
library(dplyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(pqsfinder)
library(G4SNVHunter)
library(pwalign)


# Assign genome object
hg38 <- BSgenome.Hsapiens.UCSC.hg38

# --------------------------------------------------------------------
# 2. Read your table and build GRanges
# Replace path with your file path and sheet if needed
tbl <- read_excel("testove.xlsx", sheet = 1)

# Expected columns: geneSymb (or geneSymbol), chr, strand, exon_start_plus100, exon_end_plus100
# Adjust the names if your Excel sheet has slightly different names
tbl2 <- tbl %>%
  dplyr::rename(
    gene  = geneSymbol,
    start = exon_start_plus100,
    end   = exon_end_plus100
  ) %>%
  mutate(
    chr   = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end))

tbl2 <- tbl2 %>%
  distinct(gene, .keep_all = TRUE)

print(names(tbl2))

gr <- GRanges(seqnames = tbl2$chr,
              ranges   = IRanges(start = tbl2$start, end = tbl2$end),
              strand   = tbl2$strand,
              gene     = tbl2$gene)

# --------------------------------------------------------------------
# 3. Extract sequences from hg38
seqs <- getSeq(hg38, gr)

# Handle strand: use `ifelse()` will convert class – instead use indexing
minus_idx <- which(strand(gr) == "-")
seqs_rc <- seqs  # keep class
if (length(minus_idx) > 0) {
  seqs_rc[minus_idx] <- reverseComplement(seqs[minus_idx])
}
names(seqs_rc) <- paste0(mcols(gr)$gene, "_", seqnames(gr), ":", start(gr), "-", end(gr))

# Now seqs_rc should be a DNAStringSet
class(seqs_rc)
# Then write
writeXStringSet(seqs_rc, filepath="extracted_sequences.fa", format="fasta")


# --------------------------------------------------------------------
# 4. Run pqsfinder on each sequence
# This will identify putative G4 motifs
pqs_hits <- lapply(seqs_rc, function(s) {
  pqsfinder(s)
})


#----- TO SAVE THE PQSFINDER results into a table ------

pqs_df <- do.call(bind_rows, lapply(names(pqs_hits), function(reg) {
  pvs <- pqs_hits[[reg]]
  if (length(pvs) == 0) return(NULL)
  
  # Convert to GRanges
  gr_pvs <- as(pvs, "GRanges")
  # Extract sequences for each PQS hit
  seqs_pvs <- as.character(as(pvs, "DNAStringSet"))
  
  df <- as.data.frame(gr_pvs) %>%
    mutate(
      region      = reg,
      pq_sequence = seqs_pvs
    ) %>%
    select(region, start, width, score, strand, pq_sequence)
  
  df
}))

# If you only need gene/region + sequence (and perhaps score), you can subset:
export_df <- pqs_df %>%
  select(region, pq_sequence, score)

# Write to CSV
write.csv(export_df, file = "pqsfinder_results_gene_sequence.csv", row.names = FALSE)







# 4.1 ---- G4HUNTER part:

# Given you already have:
# seqs_rc  — a DNAStringSet named by gene region
# gr       — GRanges aligned to these sequences (with gene etc)

library(G4SNVHunter)

# Define parameters — you can choose window size and threshold per Bedrat et al. 2016
window_size <- 25
threshold   <- 1.2


# Assume you have already run:
G4_detected <- G4HunterDetect(seqs_rc, threshold = threshold, window_size = window_size)

# Then export predicted G4s:
exportG4(G4 = G4_detected,
         filename = "G4hunter.csv",
         include_metadata = TRUE,
         revcomp_antisense = FALSE)

#---------------------------------------------------------------------------------------------------------------------

g4_hits <- G4HunterDetect(seqs_rc, threshold = 1.2, window_size = 25)

# Convert to data frame
df_g4 <- as.data.frame(g4_hits)

# Add region identifier (matching names(seqs_rc))
df_g4 <- df_g4 %>% 
  mutate(region = seqnames %>% as.character())

# Summarise per region
df_summary <- df_g4 %>%
  group_by(region) %>%
  summarise(
    best_score = max(max_score, na.rm=TRUE),
    best_idx   = which.max(max_score),
    best_seq   = sequence[best_idx],
    start      = min(start[best_idx]),
    end        = max(end[best_idx]),
    strand     = unique(strand[best_idx])
  ) %>% 
  ungroup()

# Export
write.csv(df_summary, file = "g4hunter_best_hits_per_region.csv", row.names = FALSE)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pqs_df <- do.call(bind_rows, lapply(names(pqs_hits), function(reg) {
  pvs <- pqs_hits[[reg]]
  if (length(pvs) == 0) return(NULL)
  # Convert to GRanges
  gr_pvs <- as(pvs, "GRanges")
  # Also get the sequences
  seqs_pvs <- as.character(as(pvs, "DNAStringSet"))
  df <- as.data.frame(gr_pvs) %>%
    mutate(region = reg,
           pq_sequence = seqs_pvs) %>%
    select(region, start, width, score, strand, pq_sequence)
  df
}))

# G4Hunter sequences
g4_df <- df_g4 %>% select(region, g4_sequence = sequence)

#------------------------------------------------------------------------------------------------

# Add an index (hit number) per region for both data sets
pqs_indexed <- pqs_df %>%
  group_by(region) %>%
  mutate(pq_hit = row_number()) %>%
  ungroup()

g4_indexed <- g4_df %>%
  group_by(region) %>%
  mutate(g4_hit = row_number()) %>%
  ungroup()

# Now perform a full join so you keep all hits from both sides
combined <- full_join(
  g4_indexed %>% select(region, g4_hit, g4_sequence),
  pqs_indexed %>% select(region, pq_hit, pq_sequence),
  by = "region"
) %>%
  arrange(region, g4_hit, pq_hit)

# Export
write.csv(combined, file = "combined_g4_pqsfinder_side_by_side.csv", row.names = FALSE)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


combined2 <- combined %>%
  rowwise() %>%
  mutate(
    pct_id = if (!is.na(g4_sequence) & !is.na(pq_sequence)) {
      aln <- pwalign::pairwiseAlignment(pattern = g4_sequence,
                                        subject = pq_sequence,
                                        type = "global",
                                        gapOpening = 5,
                                        gapExtension = 2)
      pwalign::pid(aln, type = "PID3")
    } else {
      NA_real_
    }
  ) %>%
  ungroup() %>%
  select(region, g4_hit, g4_sequence, pq_hit, pq_sequence, pct_id)

write.csv(combined2, file="combined_with_similarity.csv", row.names=FALSE)

# # *********** Enhanced similarity analysis with thresholds
# combined3 <- combined2 %>%
#   mutate(
#     match_category = case_when(
#       is.na(pct_id) ~ "No comparison",
#       pct_id >= 80 ~ "High similarity (≥80%)",
#       pct_id >= 60 ~ "Moderate similarity (60-79%)",
#       pct_id >= 40 ~ "Low similarity (40-59%)",
#       TRUE ~ "Very low similarity (<40%)"
#     )
#   )
# 
# # Summary table
# similarity_summary <- combined3 %>%
#   count(match_category) %>%
#   arrange(desc(n))
# 
# print(similarity_summary)












#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

master  <- "GGGCTTCTAGGGTCAGGTCCAGGGGGCCCGGCAGGGGTCAGGTTGG"
hit1    <- "GGTGGGCTTCTAGGGTCAGGTCCAGGGGGCCCGGCAGGG"
hit2    <- "GTGGGGTAGGGGGCCTGGCCGGGAGGATGGG"
hit3    <- "GGGTGGGCATATGGAGAGAGGGCGACCGGG"


library(pwalign)  # or Biostrings, depending on version
library(Biostrings)

aln1 <- pwalign::pairwiseAlignment(pattern = master, subject = hit1, type="global")
pattern(aln1)
subject(aln1)
compareStrings(pattern(aln1), subject(aln1))
pid(aln1, type="PID3")

aln2 <- pwalign::pairwiseAlignment(pattern = master, subject = hit2, type="global")
pattern(aln2)
subject(aln2)
compareStrings(pattern(aln2), subject(aln2))
pid(aln2, type="PID3")

aln3 <- pwalign::pairwiseAlignment(pattern = master, subject = hit3, type="global")
pattern(aln3)
subject(aln3)
compareStrings(pattern(aln3), subject(aln3))
pid(aln3, type="PID3")

#```````````````````````````````````````````````````````````



# --------------------------------------------------------------------
# Create BED file with PQS coordinates for IGV visualization
# --------------------------------------------------------------------

# Build a detailed data frame with genomic coordinates for each PQS hit
pqs_bed_df <- do.call(bind_rows, lapply(seq_along(pqs_hits), function(i) {
  pvs <- pqs_hits[[i]]
  if (length(pvs) == 0) return(NULL)
  
  # Get the corresponding GRanges entry
  gene_name <- mcols(gr)$gene[i]
  gene_chr <- as.character(seqnames(gr)[i])
  gene_start <- start(gr)[i]
  gene_end <- end(gr)[i]
  gene_strand <- as.character(strand(gr)[i])
  
  # Convert PQS to GRanges (relative positions within extracted sequence)
  gr_pvs <- as(pvs, "GRanges")
  seqs_pvs <- as.character(as(pvs, "DNAStringSet"))
  
  # Calculate absolute genomic coordinates
  # PQS positions are relative to the extracted sequence
  if (gene_strand == "-") {
    # For minus strand, coordinates are reversed
    abs_starts <- gene_end - end(gr_pvs) + 1
    abs_ends <- gene_end - start(gr_pvs) + 1
  } else {
    # For plus strand, add offset
    abs_starts <- gene_start + start(gr_pvs) - 1
    abs_ends <- gene_start + end(gr_pvs) - 1
  }
  
  # Create data frame
  data.frame(
    chr = gene_chr,
    start = abs_starts,
    end = abs_ends,
    name = paste0(gene_name, "_PQS", seq_along(pvs)),
    score = score(pvs),
    strand = gene_strand,
    gene = gene_name,
    sequence = seqs_pvs,
    stringsAsFactors = FALSE
  )
}))

# Sort by chromosome and position
pqs_bed_df <- pqs_bed_df %>%
  arrange(chr, start)

# Create standard BED format (BED6 + extra columns)
# BED format: chr start end name score strand [optional fields]
bed_output <- pqs_bed_df %>%
  mutate(
    # Ensure score is between 0-1000 for BED format
    bed_score = pmin(round(score * 10), 1000),
    # Create a descriptive name with sequence
    bed_name = paste0(name, "|", sequence, "|score:", round(score, 2))
  ) %>%
  select(chr, start, end, bed_name, bed_score, strand)

# Write main BED file (all genes combined)
write.table(bed_output, 
            file = "pqsfinder_all_genes.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("Created: pqsfinder_all_genes.bed\n")

# --------------------------------------------------------------------
# Create separate BED file for each gene
# --------------------------------------------------------------------

# Create a directory for individual gene BED files
if (!dir.exists("gene_bed_files")) {
  dir.create("gene_bed_files")
}

# Split by gene and create individual files
pqs_bed_df %>%
  group_by(gene) %>%
  group_walk(~ {
    gene_name <- .y$gene
    bed_gene <- .x %>%
      mutate(
        bed_score = pmin(round(score * 10), 1000),
        bed_name = paste0(gene_name, "_PQS", row_number(), "|", sequence, "|score:", round(score, 2))
      ) %>%
      select(chr, start, end, bed_name, bed_score, strand)
    
    filename <- file.path("gene_bed_files", paste0(gene_name, "_pqsfinder.bed"))
    write.table(bed_gene,
                file = filename,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE)
  })

cat(sprintf("Created %d individual gene BED files in 'gene_bed_files/' directory\n", 
            n_distinct(pqs_bed_df$gene)))

# --------------------------------------------------------------------
# Create a summary table with genomic coordinates
# --------------------------------------------------------------------

pqs_summary_with_coords <- pqs_bed_df %>%
  select(gene, chr, start, end, strand, score, sequence) %>%
  arrange(gene, start)

write.csv(pqs_summary_with_coords, 
          file = "pqsfinder_results_with_genomic_coords.csv",
          row.names = FALSE)

cat("Created: pqsfinder_results_with_genomic_coords.csv\n")

# --------------------------------------------------------------------
# Optional: Create BED file with color coding by score
# --------------------------------------------------------------------

# BED9 format includes color (RGB)
bed_colored <- pqs_bed_df %>%
  mutate(
    bed_score = pmin(round(score * 10), 1000),
    bed_name = paste0(name, "|", sequence),
    thickStart = start,  # BED9 format
    thickEnd = end,
    # Color gradient: low score (blue) to high score (red)
    # Normalize score to 0-1 range
    score_norm = (score - min(score)) / (max(score) - min(score)),
    itemRgb = sprintf("%d,%d,%d",
                      round(score_norm * 255),  # Red channel
                      0,                         # Green channel
                      round((1 - score_norm) * 255))  # Blue channel
  ) %>%
  select(chr, start, end, bed_name, bed_score, strand, 
         thickStart, thickEnd, itemRgb)

write.table(bed_colored,
            file = "pqsfinder_all_genes_colored.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("Created: pqsfinder_all_genes_colored.bed (with color coding by score)\n")

# --------------------------------------------------------------------
# Print summary statistics
# --------------------------------------------------------------------

cat("\n=== PQSfinder BED Export Summary ===\n")
cat(sprintf("Total PQS hits: %d\n", nrow(pqs_bed_df)))
cat(sprintf("Genes with PQS: %d\n", n_distinct(pqs_bed_df$gene)))
cat(sprintf("Average PQS per gene: %.2f\n", 
            nrow(pqs_bed_df) / n_distinct(pqs_bed_df$gene)))
cat(sprintf("Score range: %.2f - %.2f\n", 
            min(pqs_bed_df$score), max(pqs_bed_df$score)))
cat("\nFiles created:\n")
cat("  1. pqsfinder_all_genes.bed - All PQS in one file\n")
cat("  2. gene_bed_files/*.bed - Individual BED file per gene\n")
cat("  3. pqsfinder_all_genes_colored.bed - Color-coded by score\n")
cat("  4. pqsfinder_results_with_genomic_coords.csv - Detailed table\n")`

# --------------------------------------------------------------------
# Create BED file with G4Hunter coordinates for IGV visualization
# --------------------------------------------------------------------

# Assume you already have:
# gr       — your original GRanges corresponding to each region (as before)
# seqs_rc  — DNAStringSet of extracted sequences (named by region) (as before)
# g4_hits  — GRanges object returned by G4HunterDetect(seqs_rc, …) 
#             containing predicted G4s, with metadata columns including score, max_score, sequence

# Build detailed data.frame with genomic coordinates for each G4 hit
g4_bed_df <- do.call(bind_rows, lapply(seq_along(g4_hits), function(i) {
  # For each hit in g4_hits (i.e., each predicted G4 region)
  # But note: g4_hits is a single GRanges of many hits — so better split by region name
  gr_hit <- g4_hits[i]
  
  # Identify which original region this hit belongs to: easiest if seqnames(gr_hit) matches names(seqs_rc)
  region_name <- as.character(seqnames(gr_hit))
  
  # Find index in original gr of that region
  idx_orig <- match(region_name, names(seqs_rc))
  if (is.na(idx_orig)) {
    stop("Region name in g4_hits not found in original ‘seqs_rc’ names: ", region_name)
  }
  gene_name   <- mcols(gr)$gene[idx_orig]
  gene_chr    <- as.character(seqnames(gr)[idx_orig])
  gene_start  <- start(gr)[idx_orig]
  gene_end    <- end(gr)[idx_orig]
  gene_strand <- as.character(strand(gr)[idx_orig])
  
  # Coordinates relative to extracted sequence
  rel_start <- start(gr_hit)
  rel_end   <- end(gr_hit)
  
  # Compute absolute genomic coordinates
  if (gene_strand == "-") {
    abs_start <- gene_end - rel_end + 1
    abs_end   <- gene_end - rel_start + 1
  } else {
    abs_start <- gene_start + rel_start - 1
    abs_end   <- gene_start + rel_end   - 1
  }
  
  # Grab metadata columns: e.g., score, max_score, sequence
  score_val    <- mcols(gr_hit)$score
  max_score_val<- mcols(gr_hit)$max_score
  seq_val      <- mcols(gr_hit)$sequence
  
  data.frame(
    chr      = gene_chr,
    start    = abs_start,
    end      = abs_end,
    name     = paste0(gene_name, "_G4HIT", i),
    score    = score_val,
    max_score= max_score_val,
    strand   = gene_strand,
    gene     = gene_name,
    sequence = seq_val,
    stringsAsFactors = FALSE
  )
}))

# Sort by chromosome and position
g4_bed_df <- g4_bed_df %>% arrange(chr, start)

# Create standard BED format (BED6 + extra optional columns)
bed_output_g4 <- g4_bed_df %>%
  mutate(
    bed_score = pmin(round(score * 10), 1000),
    bed_name  = paste0(name, "|seq:", sequence, "|score:", round(score,2), "|max_score:", round(max_score,2))
  ) %>%
  select(chr, start, end, bed_name, bed_score, strand)

# Write main BED file (all regions combined)
write.table(bed_output_g4,
            file = "G4Hunter_all_regions.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("Created: G4Hunter_all_regions.bed\n")

# --------------------------------------------------------------------
# (Optional) Create colored BED based on max_score or score
# --------------------------------------------------------------------
bed_colored_g4 <- g4_bed_df %>%
  mutate(
    bed_score = pmin(round(score * 10), 1000),
    bed_name  = paste0(name, "|seq:", sequence),
    thickStart = start,
    thickEnd   = end,
    score_norm = (max_score - min(max_score, na.rm=TRUE)) / (max(max_score, na.rm=TRUE) - min(max_score, na.rm=TRUE)),
    itemRgb    = sprintf("%d,%d,%d",
                         round(score_norm * 255),  # red
                         0,                         # green
                         round((1 - score_norm) * 255))  # blue
  ) %>%
  select(chr, start, end, bed_name, bed_score, strand,
         thickStart, thickEnd, itemRgb)

write.table(bed_colored_g4,
            file = "G4Hunter_all_regions_colored.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("Created: G4Hunter_all_regions_colored.bed (with color coding by max_score)\n")

# --------------------------------------------------------------------
# Print summary
# --------------------------------------------------------------------
cat("\n=== G4Hunter BED Export Summary ===\n")
cat(sprintf("Total G4 hits: %d\n", nrow(g4_bed_df)))
cat(sprintf("Regions with G4 hits: %d\n", length(unique(g4_bed_df$gene))))
cat(sprintf("Average G4 hits per region: %.2f\n",
            nrow(g4_bed_df) / length(unique(g4_bed_df$gene))))
cat(sprintf("Score (score) range: %.2f - %.2f\n",
            min(g4_bed_df$score, na.rm=TRUE), max(g4_bed_df$score, na.rm=TRUE)))
cat(sprintf("Max_score range: %.2f - %.2f\n",
            min(g4_bed_df$max_score, na.rm=TRUE), max(g4_bed_df$max_score, na.rm=TRUE)))
cat("\nFiles created:\n")
cat("  1. G4Hunter_all_regions.bed\n")
cat("  2. G4Hunter_all_regions_colored.bed\n")











