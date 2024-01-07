library(readr)
library(tidyverse)
library(ggplot2)


md = read_csv("metadata.csv")
md = md %>% select(sample, population, time, clone, sample_type)

X = read_csv("all_evidence_for_plotting.csv")


#########################################################################
## Perform filtering
filtered_mutation_stats = read_csv("final_mutation_stats.csv")
ids_to_keep = filtered_mutation_stats$id

######################################
# Plot

my_theme =  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dir.create("graphs/mutation",recursive=T)
for(mut in ids_to_keep) {
  cat("Plotting mutation: ", mut, "\n")
  Y = X %>% filter(id==mut)
  
  Y$sample = factor(Y$sample, levels=Y$sample)
  Y$line = paste0("Line ", Y$line)
  Y$line[Y$line=="Line 0"] = "Ancestor"
  
  Y$line = factor(Y$line, levels = unique(Y$line))
  Y$time = factor(Y$time, levels=unique(Y$time))
  
  ggplot(Y, aes(x = sample, y = frequency, color = time, shape=sample_type)) +
    scale_color_discrete() +
    scale_x_discrete(labels=NULL) +
    scale_y_continuous(limits = c(0,1)) +
    geom_point() +
    facet_wrap(~line, drop=T, scales = "free_x") + 
    my_theme
  
  ggsave(file.path("graphs", "mutation", paste0(mut, ".pdf")), width=8, height=6)
}
