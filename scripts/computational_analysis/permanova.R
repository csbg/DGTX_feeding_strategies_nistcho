#### R script used to perform PERmutational Multivariate ANalysis of VAriance (PERMANOVA) for 
#### the DGTX feeding strategies project
# Developed by Veronika Schäpertöns
# Date created: 05/12/2026
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R version  4.3.0 (2023-04-21 ucrt) -- "Already Tomorrow"

# Libraries needed
library(tidyverse)
library(here)
library(compositions)
library(vegan)

# Data load
input_file_path <- here::here("analysis", "matrix_meta_four_br.RData")

load(file = input_file_path)
data.matrix
meta


strategy_order <- c("STD", "HiF", "LoG", "HIP", "STD+", "HIP+", "LoG+")


meta <- meta %>%
  mutate(condition_tp = paste(condition, timepoint, sep = "_")) %>%
  mutate(condition_abrev = case_when(
    condition == "A" ~ "STD",
    condition == "B" ~ "STD+",
    condition == "G" ~ "LoG",
    condition == "C" ~ "LoG+",
    condition == "D" ~ "HiF",
    condition == "E" ~ "HIP",
    condition == "F" ~ "HIP+")
  ) %>%
  mutate(phase = ifelse(timepoint == "120", "early", "late") ) %>%
  mutate(condition_abrev = factor(condition_abrev, levels = strategy_order),
         phase = factor(phase,levels = c("early", "late"))) %>%
  mutate(condition_abrev_tp = paste(condition_abrev, timepoint, sep = "_"))

meta <- meta %>%
  arrange(phase, condition_abrev)

data.matrix <- data.matrix[, meta$sample_name, drop = FALSE]

# Perform CLR transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix


# permanova ---------------------------------------------------------------

all(colnames(clr_data.matrix) == meta$sample_name)


#calculation of distance matrix (Euclidean distances due to clr-transformed data)
dist.eu <- dist(t(clr_data.matrix))

#check assumption: multivariate dispersion
#PERMANOVA assumes homogeneous group dispersions.
#If unequal, results may reflect dispersion, not centroid differences.
mod <- betadisper(dist.eu, meta$condition_abrev)
anova(mod)
plot(mod)
boxplot(mod)

mod <- betadisper(dist.eu, meta$phase)
anova(mod)
plot(mod)
boxplot(mod)


#0. univariate analysis the classical way
summary(aov(formula = clr_data.matrix ~ condition_abrev, data = meta))

#1. univariate analysis
adonis2(t(clr_data.matrix) ~ condition_abrev, data = meta, method = "euclidean")

adonis2(dist.eu ~ condition_abrev, data = meta)
adonis2(dist.eu ~ phase, data = meta)



#2. multivariate analysis
#effect of terms independently
res1 <- adonis2(dist.eu ~ phase + condition_abrev, data = meta)
res1

#terms are crossed
res2 <- adonis2(dist.eu ~ phase * condition_abrev, data = meta)
res2










