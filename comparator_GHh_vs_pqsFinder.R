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

library(pwalign)
library(readxl)
library(dplyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(pqsfinder)
library(G4SNVHunter)

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
