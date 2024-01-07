## This takes a directory of CSV files from gdtools ANNOTATE -e (with evidence)
## and splits out the info needed for RA and JC items into different output files.

library(tidyverse)

files <- list.files(path="gd_user_evidence_csv", pattern="*.csv", full.names=TRUE, recursive=FALSE)

# Ignore samples marked as having low coverage.
metadata = read_csv("metadata.csv")
dont_include_files = (metadata %>% filter(failed_low_coverage) )
length(files)
for (i in 1:nrow(dont_include_files)) {
  cat("Ignoring sample:", dont_include_files[i,]$sample, "\n")
  files = files[!grepl(dont_include_files[i,]$sample, files)]
}
length(files)

RAs = data.frame()
JCs = data.frame()

for (f in files) {
  cat("Loading:", f, "\n")
  X = read_csv(f, col_types = cols(.default = "c"))
  
  Y = X %>% filter(type %in% c("RA"))
  Y = Y %>% select(title, type, seq_id, position, insert_position, ref_base, new_base, frequency, total_cov, new_cov, ref_cov)
  Y = Y %>% mutate(id=paste(type, seq_id, position, insert_position, ref_base, new_base, sep="__"))
  RAs = RAs %>% bind_rows(Y)
  
  Y = X %>% filter(type %in% c("JC")) %>% filter(is.na(circular_chromosome)) %>% select(title, type, side_1_seq_id, side_1_position, side_1_strand, side_2_seq_id, side_2_position, side_2_strand, overlap, frequency, side_1_possible_overlap_registers, side_1_read_count, side_1_redundant, side_2_possible_overlap_registers, side_2_read_count, side_2_redundant, junction_possible_overlap_registers, new_junction_read_count)
  Y = Y %>% mutate(id=paste(type, side_1_seq_id, side_1_position, side_1_strand, side_2_seq_id, side_2_position, side_2_strand, overlap, sep="__"))
  JCs = JCs %>% bind_rows(Y)
}

## We only want to keep mutations that are in all samples...
## sometimes new mutations became visible on re-runs?

num_samples = length(unique(RAs$title))

## Only keep mutations represented in all samples
RA_counts = RAs %>% group_by(id) %>% summarize(n=n())
RA_counts = RA_counts %>% filter(n==num_samples)

RAs = RAs %>% filter(id %in% RA_counts$id) 

write_csv(RAs, "compare_RA_evidence.csv")


JC_counts = JCs %>% group_by(id) %>% summarize(n=n())
JC_counts = JC_counts %>% filter(n==num_samples)

JCs = JCs %>% filter(id %in% JC_counts$id)

write_csv(JCs, "compare_JC_evidence.csv")
