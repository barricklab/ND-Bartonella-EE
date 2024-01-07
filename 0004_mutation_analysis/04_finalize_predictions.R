library(readr)
library(tidyverse)
library(ggplot2)


X = read_csv("all_mutation_stats.csv")
cat(nrow(X))

Y = X

mutations_with_minimum_p_value = (Y %>% filter(adjusted.p.value < 1E-6))$id
cat(length(mutations_with_minimum_p_value))

mutations_with_20percent_frequency_in_population_samples = (Y %>% filter(population_frequency_max >= 0.20))$id
cat(length(mutations_with_20percent_frequency_in_population_samples))

Y$fraction_of_lines_with_samples_with_5percent_mutation = Y$lines_with_samples_with_5percent_mutation / Y$num_lines
mutations_that_are_not_in_too_many_populations = (Y %>% filter(fraction_of_lines_with_samples_with_5percent_mutation < 0.50))$id
cat(length(mutations_that_are_not_in_too_many_populations))

Y$fraction_of_samples_with_NAs = Y$samples_with_NA_frequency / Y$num_samples
mutations_that_dont_have_too_many_NA_samples = (Y %>% filter(fraction_of_lines_with_samples_with_5percent_mutation >= 0.05))$id
cat(length(mutations_that_dont_have_too_many_NA_samples))

mutations_that_have_a_high_enough_total_frequency = (Y %>% filter(sample_total_frequency >= 1))$id
cat(length(mutations_that_have_a_high_enough_total_frequency))

#final_mutations_kept_due_to_criteria = intersect(mutations_with_minimum_p_value, mutations_with_20percent_frequency_in_population_samples)

final_mutations_kept_due_to_criteria = mutations_with_minimum_p_value
cat(length(final_mutations_kept_due_to_criteria))

final_mutations_kept_due_to_criteria = intersect(final_mutations_kept_due_to_criteria, mutations_with_20percent_frequency_in_population_samples)
cat(length(final_mutations_kept_due_to_criteria))

final_mutations_kept_due_to_criteria = intersect(final_mutations_kept_due_to_criteria, mutations_that_are_not_in_too_many_populations)
cat(length(final_mutations_kept_due_to_criteria))

final_mutations_kept_due_to_criteria = intersect(final_mutations_kept_due_to_criteria, mutations_that_dont_have_too_many_NA_samples)
cat(length(final_mutations_kept_due_to_criteria))

final_mutations_kept_due_to_criteria = intersect(final_mutations_kept_due_to_criteria, mutations_that_have_a_high_enough_total_frequency)
cat(length(final_mutations_kept_due_to_criteria))

mutations_fixed_in_clonal_isolates = (Y %>% filter(clone_frequency_max >= 0.90))$id
cat(length(mutations_fixed_in_clonal_isolates))

final_mutations = union(final_mutations_kept_due_to_criteria, mutations_fixed_in_clonal_isolates)
cat(length(final_mutations))

Y = X %>% filter(id %in% final_mutations)

write_csv(Y, "final_mutation_stats.csv")


Z = read_csv("all_evidence_for_plotting.csv")
Z = Z %>% filter(id %in% final_mutations)
Z = Z %>% arrange(id, line, desc(sample_type), time, clone)
write_csv(Z, "final_mutation_frequencies.csv")



