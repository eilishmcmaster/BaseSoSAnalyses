library('ggplot2')
library('ggthemes')
library('dplyr')
library('reshape2')
library('stringr')

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/simulate_bottleneck.R?raw=TRUE")

num_loci <- 100
He_specified <- 0.32 # go slightly higher
Ho_specified <- 0.2
n <- 10

generations <- 50
pop_size <- rep(c(5,10),10) # more reps means a smoother outcome -- final data is mean 
outcross <- c(1,0.5,0) # outcross rates (0=selfing only, 1=non-self only)

out_list <- simulate_populations(pop_size, outcross, generations, num_loci, He_specified, Ho_specified)


He_df <- get_generation_df(out_list, 'He')
Ho_df <- get_generation_df(out_list, 'Ho')


ggplot(Ho_df, aes(x=generations, y=Ho, groups=factor(popsize_outcrossrate),color=factor(outcross)))+
  geom_point()+theme_bw()+geom_line()+
  geom_hline(yintercept = c(0.062, 0.038), linetype="dotted")+
  # ylim(0)+
  geom_hline(yintercept = Ho_specified, color="blue", linetype='dotted')+
  geom_hline(yintercept = He_specified, color="red", linetype='dotted')+
  facet_wrap(popsize~., drop=FALSE)+
  scale_color_brewer(palette="Set2")


ggplot(He_df, aes(x=generations, y=He, groups=factor(popsize_outcrossrate),color=factor(outcross)))+
  geom_point()+theme_bw()+geom_line()+
  geom_hline(yintercept = c(0.062, 0.038), linetype="dotted")+
  # ylim(0)+
  geom_hline(yintercept = Ho_specified, color="blue", linetype='dotted')+
  geom_hline(yintercept = He_specified, color="red", linetype='dotted')+
  facet_wrap(popsize~., drop=FALSE)+
  scale_color_brewer(palette="Set2")


