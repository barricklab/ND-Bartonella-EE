library(readr)
library(tidyverse)
library(ggplot2)

md = read_csv("metadata.csv")
md = md %>% select(sample, line, time, clone, sample_type)

JC = read_csv("compare_JC_evidence.csv")
JC = JC %>% rename(sample = title)
JC = JC %>% left_join(md, by="sample")
#Can't do anything without frequency
JC = JC %>% filter(!(side_1_redundant & side_2_redundant))

RA = read_csv("compare_RA_evidence.csv")
RA = RA %>% rename(sample = title)
RA = RA %>% left_join(md, by="sample")

#Combine RA and JC into X and sort
X = JC %>% select(sample, line, time, clone, sample_type, id, frequency)
X = X %>% bind_rows(RA %>% select(sample, line, time, clone, sample_type, id, frequency))
X = X %>% arrange(line, time, desc(sample_type), clone)

#Filter to populations for rate difference stats
JC = JC  %>% filter(sample_type=="population")
RA = RA  %>% filter(sample_type=="population")

if (!file.exists("rate_difference_tests_evidence.csv"))  {

##### Handle JC evidence

  JC_pvals = data.frame()
  for(mut in unique(JC$id)) {
    cat("Testing mutation: ", mut, "\n")
    Y = JC %>% filter(id==mut)
    
    ## Drop rows with zero total_coverage
    #side_1_redundant_multiplier = 1- Y$side_1_redundant[1]
    #side_2_redundant_multiplier = 1 - Y$side_2_redundant[1]
    
    Y$ref_junction_read_count = 0
    Y$ref_possible_overlap_registers = 0;
    
    if (!Y$side_1_redundant[1]) {
      Y$ref_junction_read_count = Y$ref_junction_read_count + Y$side_1_read_count
      Y$ref_possible_overlap_registers = Y$ref_possible_overlap_registers + Y$side_1_possible_overlap_registers
    }
    
    if (!Y$side_2_redundant[1]) {
      Y$ref_junction_read_count = Y$ref_junction_read_count + Y$side_2_read_count
      Y$ref_possible_overlap_registers = Y$ref_possible_overlap_registers + Y$side_2_possible_overlap_registers
    }
    
    Y$total_junction_read_count = Y$ref_junction_read_count + Y$new_junction_read_count
    Y$total_possible_overlap_registers = Y$ref_possible_overlap_registers + Y$junction_possible_overlap_registers
    
    ## Remove ones with no informative reads
    Y = Y %>% filter(total_junction_read_count > 0)
    
    ## The relative rates of observing new reads supporting the junction scales with a time bases
    ## proportional to the total number of reads and the different lengths of reads that could support junction
    
    Y$junction_offset = (Y$junction_possible_overlap_registers / Y$total_possible_overlap_registers) * Y$total_junction_read_count
    
    # need to expand out rows for each sample and give them 0 (ref) or 1 (variant) values
    # so we can use logistic regression
    
    fit1 = glm(new_junction_read_count~offset(log(junction_offset)), data=Y, family = poisson(link="log") )
    
    fit2 = glm(new_junction_read_count~offset(log(junction_offset))+sample, data=Y, family = poisson(link="log"), maxit=100000)
    a = anova(fit1, fit2, test="LRT")
    pval = a$"Pr(>Chi)"[2]
    
    cat("p.value = ", pval, "\n")
    
    JC_pvals = JC_pvals %>% bind_rows(data.frame(id=mut, p.value = pval))
  }  
  write_csv(JC_pvals, "JC_pvals.csv")
  
  ##### Handle RA evidence
  
  # Calculate model as to whether there are deviations from one "rate" of generating reads for the variant.
  RA = RA %>% separate(new_cov, c("new_top_cov", "new_bottom_cov"), sep="/")
  RA = RA %>% separate(total_cov, c("total_top_cov", "total_bottom_cov"), sep="/")
  RA$new_top_cov = as.numeric(RA$new_top_cov)
  RA$new_bottom_cov = as.numeric(RA$new_bottom_cov)
  RA$total_top_cov = as.numeric(RA$total_top_cov)
  RA$total_bottom_cov = as.numeric(RA$total_bottom_cov)
  RA = RA %>% mutate(new_cov = new_top_cov+new_bottom_cov, total_cov = total_top_cov+total_bottom_cov)
  RA = RA %>% select(-new_top_cov, -new_bottom_cov, -total_top_cov, -total_bottom_cov)
  
  RA_pvals = data.frame()
  for(mut in unique(RA$id)) {
    cat("Testing mutation: ", mut, "\n")
    Y = RA %>% filter(id==mut)
    
    ## Drop rows with zero total_coverage
    
    Y = Y %>% filter(total_cov > 0)
    
    Y$sample = factor(Y$sample, levels=Y$sample)
    Y$line = paste0("Line ", Y$line)
    Y$line[Y$line=="Line 0"] = "Ancestor"
    
    successes = Y$new_cov
    trials = Y$total_cov
    failures = trials - successes
    
    # need to expand out rows for each sample and give them 0 (ref) or 1 (variant) values
    # so we can use logistic regression
    
    fit1 = glm(new_cov~offset(log(total_cov)), data=Y, family = poisson(link="log") )
    
    fit2 = glm(new_cov~offset(log(total_cov))+sample, data=Y, family = poisson(link="log"), maxit=100000)
    a = anova(fit1, fit2, test="LRT")
    pval = a$"Pr(>Chi)"[2]
    
    cat("p.value = ", pval, "\n")
    
    RA_pvals = RA_pvals %>% bind_rows(data.frame(id=mut, p.value = pval))
  }  
  write_csv(RA_pvals, "RA_pvals.csv")
  
  #### adjust p-values across all RA and JC samples
  
  all_pvals = JC_pvals %>% bind_rows(RA_pvals)
  all_pvals$adjusted.p.value = p.adjust(all_pvals$p.value, method="BH")
  
  write_csv(all_pvals, "rate_difference_tests_evidence.csv")
}
