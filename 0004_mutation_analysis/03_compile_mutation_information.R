library(readr)
library(tidyverse)
library(ggplot2)

md = read_csv("metadata.csv")
md = md %>% select(sample, line, time, clone, sample_type)

### Create file with all data merged together if we haven't already
if (!file.exists("all_evidence_for_plotting.csv")) {
  
  JC = read_csv("compare_JC_evidence.csv")
  JC = JC %>% rename(sample = title)
  JC = JC %>% left_join(md, by="sample")
  #Can't do anything without frequency
  JC = JC %>% filter(!(side_1_redundant & side_2_redundant))
  JC = JC %>% mutate(total_reads = side_1_read_count + side_2_read_count + new_junction_read_count)
  
  RA = read_csv("compare_RA_evidence.csv")
  RA = RA %>% rename(sample = title)
  RA = RA %>% left_join(md, by="sample")
  RA = RA %>% separate(total_cov, c("total_top_cov", "total_bottom_cov"), sep="/")
  RA$total_top_cov = as.numeric(RA$total_top_cov)
  RA$total_bottom_cov = as.numeric(RA$total_bottom_cov)
  RA = RA %>% mutate(total_reads = total_top_cov+total_bottom_cov)
  
  #Combine RA and JC into X and sort
  X = JC %>% select(sample, line, time, clone, sample_type, id, frequency, total_reads)
  X = X %>% bind_rows(RA %>% select(sample, line, time, clone, sample_type, id, frequency, total_reads))
  X = X %>% arrange(line, time, desc(sample_type), clone)
  
  X$frequency[X$total_reads<10] = NA
  X$frequency[is.na(X$total_reads)] = NA
  
  write_csv(X, "all_evidence_for_plotting.csv")
  
} else {
  X = read_csv("all_evidence_for_plotting.csv")
}

# General dataset stats....
num_lines = length(unique((X %>% filter(line!=0))$line))
num_samples = length(unique((X %>% filter(line!=0))$sample))
num_clones = length(unique((X %>% filter(line!=0,sample_type=="clone"))$sample))
num_populations = length(unique((X %>% filter(line!=0,sample_type=="population"))$sample))

###############################################################################
#### Reload stats generated in previous steps

mutation_stats = data.frame(id=unique(X$id))

all_pvals = read_csv("rate_difference_tests_evidence.csv")
mutation_stats = mutation_stats %>% left_join(all_pvals, by="id")

#all_sad = read_csv("time_course_difference_test.csv")
#mutation_stats = mutation_stats %>% left_join(all_sad, by="id")

write_csv(mutation_stats, "all_mutation_stats.csv")

###############################################################################
#### Some junctions may have lots of NA freqs b/c there are no reads
#### All of these are included in failing the population Poisson regression test

samples_with_NA_frequency = X %>% 
  filter(line!=0) %>% 
  group_by(id) %>%
  filter(is.na(frequency)) %>% 
  summarize(samples_with_NA_frequency=n())

mutation_stats = mutation_stats %>% left_join(samples_with_NA_frequency, by="id")
mutation_stats$samples_with_NA_frequency[is.na(mutation_stats$samples_with_NA_frequency)] = 0

###############################################################################
#### High frequency in clonal samples

clones_with_100percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id) %>%
  filter(frequency==1) %>% 
  summarize(clones_with_100percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(clones_with_100percent_mutation, by="id")
mutation_stats$clones_with_100percent_mutation[is.na(mutation_stats$clones_with_100percent_mutation)] = 0

clones_with_90percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id) %>%
  filter(frequency>=0.9) %>% 
  summarize(clones_with_90percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(clones_with_90percent_mutation, by="id")
mutation_stats$clones_with_90percent_mutation[is.na(mutation_stats$clones_with_90percent_mutation)] = 0

clones_with_80percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id) %>%
  filter(frequency>=0.8) %>% 
  summarize(clones_with_80percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(clones_with_80percent_mutation, by="id")
mutation_stats$clones_with_80percent_mutation[is.na(mutation_stats$clones_with_80percent_mutation)] = 0

lines_with_clones_with_100percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id, line) %>%
  filter(frequency==1) %>% 
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_clones_with_100percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_clones_with_100percent_mutation, by="id")
mutation_stats$lines_with_clones_with_100percent_mutation[is.na(mutation_stats$lines_with_clones_with_100percent_mutation)] = 0

lines_with_clones_with_90percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id, line) %>%
  filter(frequency>=0.9) %>% 
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_clones_with_90percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_clones_with_90percent_mutation, by="id")
mutation_stats$lines_with_clones_with_90percent_mutation[is.na(mutation_stats$lines_with_clones_with_90percent_mutation)] = 0

lines_with_clones_with_80percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id, line) %>%
  filter(frequency>=0.8) %>% 
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_clones_with_80percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_clones_with_80percent_mutation, by="id")
mutation_stats$lines_with_clones_with_80percent_mutation[is.na(mutation_stats$lines_with_clones_with_80percent_mutation)] = 0

clones_with_10to90percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id) %>%
  filter(frequency>0.1, frequency<0.9) %>% 
  summarize(clones_with_10to90percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(clones_with_10to90percent_mutation, by="id")
mutation_stats$clones_with_10to90percent_mutation[is.na(mutation_stats$clones_with_10to90percent_mutation)] = 0

clones_with_20to80percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id) %>%
  filter(frequency>0.2, frequency<0.8) %>% 
  summarize(clones_with_20to80percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(clones_with_20to80percent_mutation, by="id")
mutation_stats$clones_with_20to80percent_mutation[is.na(mutation_stats$clones_with_20to80percent_mutation)] = 0

lines_with_clones_with_10to90percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id, line) %>%
  filter(frequency>0.1, frequency<0.9) %>% 
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_clones_with_10to90percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_clones_with_10to90percent_mutation, by="id")
mutation_stats$lines_with_clones_with_10to90percent_mutation[is.na(mutation_stats$lines_with_clones_with_10to90percent_mutation)] = 0

lines_with_clones_with_20to80percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id, line) %>%
  filter(frequency>0.2, frequency<0.8) %>% 
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_clones_with_20to80percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_clones_with_20to80percent_mutation, by="id")
mutation_stats$lines_with_clones_with_20to80percent_mutation[is.na(mutation_stats$lines_with_clones_with_20to80percent_mutation)] = 0

write_csv(mutation_stats, "all_mutation_stats.csv")


#######################################################################
### What is the frequency range, min, and max of the mutation?

population_mutation_frequency_info = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="population") %>% 
  filter(!is.na(frequency)) %>% 
  group_by(id) %>%
  summarize(
    population_total_frequency = sum(frequency),
    population_frequency_min = min(frequency, na.rm=TRUE),
    population_frequency_max = max(frequency, na.rm=TRUE),
    population_frequency_range = population_frequency_max - population_frequency_min
  ) 

mutation_stats = mutation_stats %>% left_join(population_mutation_frequency_info, by="id")
mutation_stats$population_total_frequency[is.na(mutation_stats$population_total_frequency)] = 0
mutation_stats$population_frequency_min[is.na(mutation_stats$population_frequency_min)] = 0
mutation_stats$population_frequency_max[is.na(mutation_stats$population_frequency_max)] = 0
mutation_stats$population_frequency_range[is.na(mutation_stats$population_frequency_range)] = 0

clone_mutation_frequency_info = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  filter(!is.na(frequency)) %>% 
  group_by(id) %>%
  summarize(
    clone_total_frequency = sum(frequency),
    clone_frequency_min = min(frequency, na.rm=TRUE),
    clone_frequency_max = max(frequency, na.rm=TRUE),
    clone_frequency_range = clone_frequency_max - clone_frequency_min
  ) 

mutation_total_frequency = X %>% 
  group_by(id) %>% 
  summarize(total_frequency=sum(frequency))

mutation_stats = mutation_stats %>% left_join(clone_mutation_frequency_info, by="id")
mutation_stats$clone_total_frequency[is.na(mutation_stats$clone_total_frequency)] = 0
mutation_stats$clone_frequency_min[is.na(mutation_stats$clone_frequency_min)] = 0
mutation_stats$clone_frequency_max[is.na(mutation_stats$clone_frequency_max)] = 0
mutation_stats$clone_frequency_range[is.na(mutation_stats$clone_frequency_range)] = 0

sample_mutation_frequency_info = X %>% 
  filter(line!=0) %>% 
  filter(!is.na(frequency)) %>% 
  group_by(id) %>%
  summarize(
    sample_total_frequency = sum(frequency),
    sample_frequency_min = min(frequency, na.rm=TRUE),
    sample_frequency_max = max(frequency, na.rm=TRUE),
    sample_frequency_range = sample_frequency_max - sample_frequency_min
  ) 

mutation_stats = mutation_stats %>% left_join(sample_mutation_frequency_info, by="id")
mutation_stats$sample_total_frequency[is.na(mutation_stats$sample_total_frequency)] = 0
mutation_stats$sample_frequency_min[is.na(mutation_stats$sample_frequency_min)] = 0
mutation_stats$sample_frequency_max[is.na(mutation_stats$sample_frequency_max)] = 0
mutation_stats$sample_frequency_range[is.na(mutation_stats$sample_frequency_range)] = 0

write_csv(mutation_stats, "all_mutation_stats.csv")

#######################################################################
### How many different population samples and lines had the mutation?

populations_with_5percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="population") %>% 
  group_by(id) %>%
  filter(frequency>=0.05) %>% 
  summarize(populations_with_5percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(populations_with_5percent_mutation, by="id")
mutation_stats$populations_with_5percent_mutation[is.na(mutation_stats$populations_with_5percent_mutation)] = 0

lines_with_populations_with_5percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="population") %>% 
  filter(frequency>=0.05) %>% 
  group_by(id, line) %>%
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_populations_with_5percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_populations_with_5percent_mutation, by="id")
mutation_stats$lines_with_populations_with_5percent_mutation[is.na(mutation_stats$lines_with_populations_with_5percent_mutation)] = 0

populations_with_10percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="population") %>% 
  group_by(id) %>%
  filter(frequency>=0.10) %>% 
  summarize(populations_with_10percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(populations_with_10percent_mutation, by="id")
mutation_stats$populations_with_10percent_mutation[is.na(mutation_stats$populations_with_10percent_mutation)] = 0

lines_with_populations_with_10percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="population") %>% 
  filter(frequency>=0.1) %>% 
  group_by(id, line) %>%
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_populations_with_10percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_populations_with_10percent_mutation, by="id")
mutation_stats$lines_with_populations_with_10percent_mutation[is.na(mutation_stats$lines_with_populations_with_10percent_mutation)] = 0

samples_with_5percent_mutation = X %>% 
  filter(line!=0) %>% 
  group_by(id) %>%
  filter(frequency>=0.05) %>% 
  summarize(samples_with_5percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(samples_with_5percent_mutation, by="id")
mutation_stats$samples_with_5percent_mutation[is.na(mutation_stats$samples_with_5percent_mutation)] = 0

lines_with_samples_with_5percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(frequency>=0.05) %>% 
  group_by(id, line) %>%
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_samples_with_5percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_samples_with_5percent_mutation, by="id")
mutation_stats$lines_with_samples_with_5percent_mutation[is.na(mutation_stats$lines_with_samples_with_5percent_mutation)] = 0

write_csv(mutation_stats, "all_mutation_stats.csv")


#######################################################################
### How many different clonal samples and line had the mutation?

clones_with_5percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id, line) %>%
  filter(frequency>=0.05) %>% 
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(clones_with_5percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(clones_with_5percent_mutation, by="id")
mutation_stats$clones_with_5percent_mutation[is.na(mutation_stats$clones_with_5percent_mutation)] = 0

lines_with_clones_with_5percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  filter(frequency>=0.05) %>% 
  group_by(id, line) %>%
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_clones_with_5percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_clones_with_5percent_mutation, by="id")
mutation_stats$lines_with_clones_with_5percent_mutation[is.na(mutation_stats$lines_with_clones_with_5percent_mutation)] = 0

clones_with_10percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  group_by(id, line) %>%
  filter(frequency>=0.1) %>% 
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(clones_with_10percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(clones_with_10percent_mutation, by="id")
mutation_stats$clones_with_10percent_mutation[is.na(mutation_stats$clones_with_10percent_mutation)] = 0

lines_with_clones_with_10percent_mutation = X %>% 
  filter(line!=0) %>% 
  filter(sample_type=="clone") %>% 
  filter(frequency>=0.1) %>% 
  group_by(id, line) %>%
  summarize(n = n()) %>%
  group_by(id) %>%
  summarize(lines_with_clones_with_10percent_mutation=n())

mutation_stats = mutation_stats %>% left_join(lines_with_clones_with_10percent_mutation, by="id")
mutation_stats$lines_with_clones_with_10percent_mutation[is.na(mutation_stats$lines_with_clones_with_10percent_mutation)] = 0


mutation_stats = mutation_stats %>% 
  mutate(
    num_lines = num_lines,
    num_samples = num_samples, 
    num_clones = num_clones,
    num_populations = num_populations
  )

write_csv(mutation_stats, "all_mutation_stats.csv")

