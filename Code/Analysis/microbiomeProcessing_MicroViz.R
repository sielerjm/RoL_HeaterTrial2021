
# Analysis Helper ---------------------------------------------------------
# Description: Instantiates important variables and settings for analyses


# Naming Conventions ------------------------------------------------------
# Description: format for naming variables/functions consistently
#   - If multiple ords per `<>`, Capitalize subsequent words (eg. <moreThanOneWord>)

# VARIABLES: 
#   <highLevelVarType>.<subtype>.<subset>.<variables>.<Modifications>
#
#   dt.<subtype>.<subset>.<modifications>  # datatables
#   df.<subtype>.<subset>.<modifications> # dataframes
#   plot.<subset>.<y-var>.<x-vars...>  # plots/figures
#   table.<subset>.<y-var>.<x-vars...>  # tables
#   mod.<subtype>.<subset>.<y>.<x-vars...>  # models (lm, glm, etc.)
#   ps.<subtype>.<subset>  # phyloseq objects

# FUNCTIONS:
#   Should be descriptive, but short enough to know the main task

# Set Environmental Variables ---------------------------------------------

# Analysis ID
analysis.ID <- paste0(
  "RoL-HeaterTrial-v2_",  # Data subsetted? If so, how? "<name>_"
  Sys.Date(),  # Date of analysis
  "_rf"
)

# Number of cores
# - Automated
#   - detectCores()  # Tells you how many cores your comp has available
#   - If cores greater than 4, use ~90%. If 4 or less, use ~50%
num.cores = ifelse(detectCores() > 4, round(detectCores()*.9), round(detectCores()*.5) )
# num.cores = 1


# Import Data -------------------------------------------------------------

# Load Phyloseq Object
# ps.all <- readRDS(paste0(path.input, "/phyloseq_rarefied_2022-08-17.rds"))
# ps.all <- ps.rar

# ps.all <- readRDS("~/Dropbox/Mac (2)/Documents/Sharpton_Lab/Projects_Repository/Rules_of_Life/RoL_HeaterTrial2021/Data/Output/dada2_1.24.0_2022-09-12/phyloseq_cleaned_rarefied_2022-09-12.rds")

ps.all <- ps.rar

# sample_data(ps.all) <- sample.data.frame(ps.all) %>%
#   rename(Body.Condition.Score = Body.Condition.Trad)


# Clean data --------------------------------------------------------------

# sample_data(ps.all) <- sample.data.frame(ps.all) %>%
#     mutate(Temperature = as.factor(Temperature),
#            DPE = as.factor(DPE))



# Load Sample Data
df.all <- sample.data.frame(ps.all)
dt.all <- sample.data.table(ps.all) 

# Remove Duplicate column names
dt.all <- dt.all %>%
  select(-1) # Sometimes sample column gets duplicated going from phyloseq object to a datatable object


# Subset Data -------------------------------------------------------------
# - Different analyses will required looking at different subsets of the data
#   - All at T0: Technically all controls. 
#       - How do diets differ?
#   - Controls, T0-T1: 
#       - How do diets differ across development?
#   - T1: 
#       - How do diets differ across treatments?
#   - Exposed/Final: 
#       - How does Diet impact X, depending on Treatment/Infection?
#       - How does Infection.Status impact X?


## Initial ----------------------------------------------------------

df.T0 <- df.all %>%
  filter(DPE == 0)

# Datatable
dt.T0 <- dt.all %>%
  filter(DPE == 0)

# Phyloseq Object
ps.T0 <- ps.all %>%
  subset_samples(DPE == 0)

## Control: Initial to Final ----------------------------------------------------

# Dataframe
df.conT0TF <- df.all %>%
  filter(Treatment == "Control" | (Treatment == "Exposed" & DPE == 0))

# Datatable
dt.conT0TF <- dt.all %>%
  filter(Treatment == "Control"| (Treatment == "Exposed" & DPE == 0))

# Phyloseq Object
ps.conT0TF <- ps.all %>%
  subset_samples(Treatment == "Control"| (Treatment == "Exposed" & DPE == 0))

## Final ----------------------------------------------------------------
# Dataframe
df.TF <- df.all %>%
  filter(DPE == 42)

# Datatable
dt.TF <- dt.all %>%
  filter(DPE == 42)

# Phyloseq Object
ps.TF <- ps.all %>%
  subset_samples(DPE == 42)

## Final and Control ----------------------------------------------------------------
# Dataframe
df.conTF <- df.all %>%
  filter(DPE == 42 & Treatment == "Control")

# Datatable
dt.conTF <- dt.all %>%
  filter(DPE == 42 & Treatment == "Control")

# Phyloseq Object
ps.conTF <- ps.all %>%
  subset_samples(DPE == 42 & Treatment == "Control")

## Final and Exposed ------------------------------------------------
# - See what the effect of exposure had on final time point 

# Dataframe
df.expTF <- df.all %>%
  filter(DPE == 42 & Treatment == "Exposed")

# Datatable
dt.expTF <- dt.all %>%
  filter(DPE == 42 & Treatment == "Exposed")

# Phyloseq Object
ps.expTF <- ps.all %>%
  subset_samples(DPE == 42 & Treatment == "Exposed")

## Exposed: Initial to Final ----------------------------------------------------

# Dataframe
df.expT0TF <- df.all %>%
  filter(Treatment == "Exposed")

# Datatable
dt.expT0TF <- dt.all %>%
  filter(Treatment == "Exposed")

# Phyloseq Object
ps.expT0TF <- ps.all %>%
  subset_samples(Treatment == "Exposed")

## Temperature: 28  ----------------------------------------------------

# Dataframe
df.temp28 <- df.all %>%
  filter(Temperature == 28)

# Datatable
dt.temp28 <- dt.all %>%
  filter(Temperature == 28)

# Phyloseq Object
ps.temp28 <- ps.all %>%
  subset_samples(Temperature == 28)

## Temperature: 32  ----------------------------------------------------

# Dataframe
df.temp32 <- df.all %>%
  filter(Temperature == 32)

# Datatable
dt.temp32 <- dt.all %>%
  filter(Temperature == 32)

# Phyloseq Object
ps.temp32 <- ps.all %>%
  subset_samples(Temperature == 32)

## Temperature: 35  ----------------------------------------------------

# Dataframe
df.temp35 <- df.all %>%
  filter(Temperature == 35)

# Datatable
dt.temp35 <- dt.all %>%
  filter(Temperature == 35)

# Phyloseq Object
ps.temp35 <- ps.all %>%
  subset_samples(Temperature == 35)

# Alpha-Diversity ---------------------------------------------------------







# Beta-Diversity ----------------------------------------------------------







# Differential Abundance --------------------------------------------------




















# [OLD] -------------------------------------------------------------------



methods.alpha <- c("Observed", "Shannon", "Simpson", "Phylogenetic") %>% 
  purrr::set_names()


## Calculate Alpha Scores -------------------------------------------------

# Creates a datatable of alpha diversity scores for each sample

dt.alphaScores.all <- alpha_base(ps.all,  # Phyloseq object
                                 methods.alpha,  # list of alpha methods
                                 "Sample",  # Column name for sample IDs
                                 T  # Set T if you have phylogenetic data
                                 ) 

                                 
dt.alphaScores.T0 <- alpha_base(ps.T0,  # Phyloseq object
                                methods.alpha,  # list of alpha methods
                                "Sample",  # Column name for sample IDs
                                T  # Set T if you have phylogenetic data
                                ) 

dt.alphaScores.conT0TF <- alpha_base(ps.conT0TF,  # Phyloseq object
                                     methods.alpha,  # list of alpha methods
                                     "Sample",  # Column name for sample IDs
                                     T  # Set T if you have phylogenetic data
) 

dt.alphaScores.TF <- alpha_base(ps.TF,  # Phyloseq object
                                methods.alpha,  # list of alpha methods
                                "Sample",  # Column name for sample IDs
                                T  # Set T if you have phylogenetic data
) 

dt.alphaScores.conTF <- alpha_base(ps.conTF,  # Phyloseq object
                                   methods.alpha,  # list of alpha methods
                                   "Sample",  # Column name for sample IDs
                                   T  # Set T if you have phylogenetic data
) 

dt.alphaScores.expTF <- alpha_base(ps.expTF,  # Phyloseq object
                                    methods.alpha,  # list of alpha methods
                                    "Sample",  # Column name for sample IDs
                                    T  # Set T if you have phylogenetic data
) 

dt.alphaScores.expT0TF <- alpha_base(ps.expT0TF,  # Phyloseq object
                                   methods.alpha,  # list of alpha methods
                                   "Sample",  # Column name for sample IDs
                                   T  # Set T if you have phylogenetic data
) 

dt.alphaScores.expT0TF <- alpha_base(ps.expT0TF,  # Phyloseq object
                                     methods.alpha,  # list of alpha methods
                                     "Sample",  # Column name for sample IDs
                                     T  # Set T if you have phylogenetic data
) 

dt.alphaScores.temp28 <- alpha_base(ps.temp28,  # Phyloseq object
                                    methods.alpha,  # list of alpha methods
                                    "Sample",  # Column name for sample IDs
                                    T  # Set T if you have phylogenetic data
) 

dt.alphaScores.temp32 <- alpha_base(ps.temp32,  # Phyloseq object
                                    methods.alpha,  # list of alpha methods
                                    "Sample",  # Column name for sample IDs
                                    T  # Set T if you have phylogenetic data
) 

dt.alphaScores.temp35 <- alpha_base(ps.temp35,  # Phyloseq object
                                     methods.alpha,  # list of alpha methods
                                     "Sample",  # Column name for sample IDs
                                     T  # Set T if you have phylogenetic data
) 

## Normalize Alpha Scores -------------------------------------------------

# Normalize alpha scores from 0 to 1
dt.alphaScores.norm.all <- norm_alpha_score(dt.alphaScores.all, df.all, methods.alpha)
dt.alphaScores.norm.T0 <- norm_alpha_score(dt.alphaScores.T0, df.T0, methods.alpha)
dt.alphaScores.norm.conT0TF <- norm_alpha_score(dt.alphaScores.conT0TF, df.conT0TF, methods.alpha)
dt.alphaScores.norm.TF <- norm_alpha_score(dt.alphaScores.TF, df.TF, methods.alpha)
dt.alphaScores.norm.conTF <- norm_alpha_score(dt.alphaScores.conTF, df.conTF, methods.alpha)
dt.alphaScores.norm.expTF <- norm_alpha_score(dt.alphaScores.expTF, df.expTF, methods.alpha)
dt.alphaScores.norm.expT0TF <- norm_alpha_score(dt.alphaScores.expT0TF, df.expT0TF, methods.alpha)
dt.alphaScores.norm.temp28 <- norm_alpha_score(dt.alphaScores.temp28, df.temp28, methods.alpha)
dt.alphaScores.norm.temp32 <- norm_alpha_score(dt.alphaScores.temp32, df.temp32, methods.alpha)
dt.alphaScores.norm.temp35 <- norm_alpha_score(dt.alphaScores.temp35, df.temp35, methods.alpha)


## Alpha Datatable -------------------------------------------------------


# Make a datatabe containing sample data and alpha diversity
dt.alphaPlus.all <- alpha_dataTable(dt.all, dt.alphaScores.norm.all)
dt.alphaPlus.T0 <- alpha_dataTable(dt.T0, dt.alphaScores.norm.T0)
dt.alphaPlus.conT0TF <- alpha_dataTable(dt.conT0TF, dt.alphaScores.norm.conT0TF)
dt.alphaPlus.TF <- alpha_dataTable(dt.TF, dt.alphaScores.norm.TF)
dt.alphaPlus.conTF <- alpha_dataTable(dt.conTF, dt.alphaScores.norm.conTF)
dt.alphaPlus.expTF <- alpha_dataTable(dt.expTF, dt.alphaScores.norm.expTF)
dt.alphaPlus.expT0TF <- alpha_dataTable(dt.expT0TF, dt.alphaScores.norm.expT0TF)
dt.alphaPlus.temp28 <- alpha_dataTable(dt.temp28, dt.alphaScores.norm.temp28)
dt.alphaPlus.temp32 <- alpha_dataTable(dt.temp32, dt.alphaScores.norm.temp32)
dt.alphaPlus.temp35 <- alpha_dataTable(dt.temp35, dt.alphaScores.norm.temp35)


# Melt (pivot_longer) data table for easy plotting and statistical analysis
dt.alphaPlus.all.melt <- melt_alphaDataTable(dt.alphaPlus.all)
dt.alphaPlus.T0.melt <- melt_alphaDataTable(dt.alphaPlus.T0)
dt.alphaPlus.conT0TF.melt <- melt_alphaDataTable(dt.alphaPlus.conT0TF)
dt.alphaPlus.TF.melt <- melt_alphaDataTable(dt.alphaPlus.TF)
dt.alphaPlus.conTF.melt <- melt_alphaDataTable(dt.alphaPlus.conTF)
dt.alphaPlus.expT0TF.melt <- melt_alphaDataTable(dt.alphaPlus.expT0TF)
dt.alphaPlus.temp28.melt <- melt_alphaDataTable(dt.alphaPlus.temp28)
dt.alphaPlus.temp32.melt <- melt_alphaDataTable(dt.alphaPlus.temp32)
dt.alphaPlus.temp35.melt <- melt_alphaDataTable(dt.alphaPlus.temp35)


# Beta-Diversity ----------------------------------------------------------

# Distance lists
distList.all <- gen.dist.matrices(ps.all, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
methods.beta <- names(distList.all) %>% set_names(., .)

distList.T0 <- gen.dist.matrices(ps.T0, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.conT0TF <- gen.dist.matrices(ps.conT0TF, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.TF <- gen.dist.matrices(ps.TF, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.conTF <- gen.dist.matrices(ps.conTF, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.expTF <- gen.dist.matrices(ps.expTF, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.expT0TF <- gen.dist.matrices(ps.expT0TF, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.temp28 <- gen.dist.matrices(ps.temp28, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.temp32 <- gen.dist.matrices(ps.temp32, methods = c("taxonomic", "phylogenetic"), cores = num.cores)
distList.temp35 <- gen.dist.matrices(ps.temp35, methods = c("taxonomic", "phylogenetic"), cores = num.cores)

## Additional
# distList.conT0T1.Gemma <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "Gemma"), methods = "taxonomic", cores = num.cores)
# distList.conT0T1.Watts <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "Watts"), methods = "taxonomic", cores = num.cores)
# distList.conT0T1.ZIRC <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "ZIRC"), methods = "taxonomic", cores = num.cores)
# distList.conT1.ZIRC <- gen.dist.matrices(subset_samples(ps.conT1, Diet == "ZIRC"), methods = "taxonomic", cores = num.cores)

# Differential Abundance ------------------------------------------------

# 
# ## Diet --------------------------------------------------------------------
# pseq = ps.conT0T1
# data = dt.conT0T1
# 
# # Relevel
# 
# 
# levels(data$Timepoint) <- factor(levels(data$Timepoint), levels = c("3pmf", "6mpf"))
# data$Timepoint <- relevel(factor(data$Timepoint), ref = "3mpf")
# 
# levels(sample_data(pseq)$Timepoint) <- factor(levels(sample_data(pseq)$Timepoint),levels = c("3pmf", "6mpf"))
# sample_data(pseq)$Timepoint <- relevel(factor(sample_data(pseq)$Timepoint), ref = "3mpf")
# 
# # ADD INTERACTION
# sample_data(pseq) <- microbiome::meta(pseq) %>%
#   mutate(Diet.Time = paste0(Diet,".", Timepoint), .after = Age) 
# 
# data <- data %>%
#   mutate(Diet.Time = paste0(Diet,".", Timepoint), .after = Age)
# 
# # Sanity check counts are correct
# 
# microbiome::meta(pseq) %>%
#   group_by(Diet.Time) %>%
#   count()
# 
# 
# # Genus level data
# Genus_data = microbiome::aggregate_taxa(pseq, "Genus")
# 
# 
# tax.level = Genus_data #phylum_data
# tax.label = "Genus"
# tax.data = taxa.data.table(tax.level)
# 
# 
# out = ancombc(phyloseq = tax.level, 
#               formula = "Diet*Timepoint", 
#               p_adj_method = "BH", 
#               # lib_cut = 0, 
#               # prv_cut = 0,
#               group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#               # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#               neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#               tol = 1e-05, 
#               max_iter = 100,
#               conserve = T, 
#               alpha = 0.05, 
#               global = T
#               #assay_name = "counts"
# )
# 
# res.TIME = out$res
# res_global.TIME = out$res_global
# 
# 
# 
# # Global
# 
# tab_w = res_global.TIME[, "W", drop = FALSE]
# tab_p = res_global.TIME[, "p_val", drop = FALSE]
# tab_q = res_global.TIME[, "q_val", drop = FALSE]
# tab_diff = res_global.TIME[, "diff_abn", drop = FALSE]
# 
# 
# # Log transform
# 
# samp_frac = out$samp_frac
# # Replace NA with 0
# samp_frac[is.na(samp_frac)] = 0 
# 
# # Add pesudo-count (1) to avoid taking the log of 0
# log_obs_abn = log(abundances(tax.level) + 1) 
# 
# # Adjust the log observed abundances
# log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
# 
# # Prep for visualization
# 
# sig_taxa = res_global %>%
#   tibble::rownames_to_column("Taxon") %>%
#   dplyr::filter(diff_abn == TRUE) %>%
#   .$Taxon
# 
# df_sig.TIME = as.data.frame(t(log_obs_abn_adj[sig_taxa, ])) %>%
#   tibble::rownames_to_column("Sample") %>%
#   dplyr::left_join(data %>%
#                      dplyr::select(Sample, Diet, Timepoint, Diet.Time),
#                    by = "Sample") %>%
#   dplyr::filter(!is.na(Diet)) %>%
#   tidyr::pivot_longer(cols = -one_of("Sample", "Diet", "Timepoint", "Diet.Time"), 
#                       names_to = "Taxon", values_to = "value")
# 
# 
# 
# 
# ## Exposure ----------------------------------------------------------------
# 
# pseq = ps.all
# data = dt.all
# 
# # Relevel
# 
# 
# levels(data$PrePostExp) <- factor(levels(data$PrePostExp), levels = c("Pre-exposure", "Unexposed", "Exposed"))
# data$PrePostExp <- relevel(factor(data$PrePostExp), ref = "Unexposed")
# 
# levels(sample_data(pseq)$PrePostExp) <- factor(levels(sample_data(pseq)$PrePostExp),levels = c("Pre-exposure", "Unexposed", "Exposed"))
# sample_data(pseq)$PrePostExp <- relevel(factor(sample_data(pseq)$PrePostExp), ref = "Unexposed")
# 
# # ADD INTERACTION
# sample_data(pseq) <- microbiome::meta(pseq) %>%
#   mutate(Diet.Exp = paste0(Diet,".", PrePostExp), .after = Age) 
# 
# data <- data %>%
#   mutate(Diet.Exp = paste0(Diet,".", PrePostExp), .after = Age)
# 
# # Sanity check counts are correct
# 
# microbiome::meta(pseq) %>%
#   group_by(Diet.Exp) %>%
#   count()
# 
# 
# # Genus level data
# Genus_data = microbiome::aggregate_taxa(pseq, "Genus")
# 
# 
# tax.level = Genus_data #phylum_data
# tax.label = "Genus"
# tax.data = taxa.data.table(tax.level)
# 
# 
# out = ancombc(phyloseq = tax.level, 
#               formula = "Diet*PrePostExp", 
#               p_adj_method = "BH", 
#               # lib_cut = 0, 
#               # prv_cut = 0,
#               group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#               # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#               neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#               tol = 1e-05, 
#               max_iter = 100,
#               conserve = T, 
#               alpha = 0.05, 
#               global = T
#               #assay_name = "counts"
#               )
# 
# res.EXP = out$res
# res_global.EXP = out$res_global
# 
# 
# 
# # Global
# 
# tab_w = res_global.EXP[, "W", drop = FALSE]
# tab_p = res_global.EXP[, "p_val", drop = FALSE]
# tab_q = res_global.EXP[, "q_val", drop = FALSE]
# tab_diff = res_global.EXP[, "diff_abn", drop = FALSE]
# 
# 
# # Log transform
# 
# samp_frac = out$samp_frac
# # Replace NA with 0
# samp_frac[is.na(samp_frac)] = 0 
# 
# # Add pesudo-count (1) to avoid taking the log of 0
# log_obs_abn = log(abundances(tax.level) + 1) 
# 
# # Adjust the log observed abundances
# log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
# 
# # Prep for visualization
# 
# sig_taxa = res_global %>%
#   tibble::rownames_to_column("Taxon") %>%
#   dplyr::filter(diff_abn == TRUE) %>%
#   .$Taxon
# 
# df_sig.EXP = as.data.frame(t(log_obs_abn_adj[sig_taxa, ])) %>%
#   tibble::rownames_to_column("Sample") %>%
#   dplyr::left_join(data %>%
#                      dplyr::select(Sample, Diet, PrePostExp, Diet.Exp),
#                    by = "Sample") %>%
#   dplyr::filter(!is.na(Diet)) %>%
#   tidyr::pivot_longer(cols = -one_of("Sample", "Diet", "PrePostExp", "Diet.Exp"), 
#                       names_to = "Taxon", values_to = "value")
# 
# 
# 
# 
# 
# ## Diet:Time -----------------------------------------------------------
# 
# # Load Data
# 
# pseq.Gemma = ps.conT0T1 %>% subset_samples(Diet == "Gemma")
# data.Gemma = dt.conT0T1 %>% filter(Diet == "Gemma")
# 
# pseq.Watts = ps.conT0T1 %>% subset_samples(Diet == "Watts")
# data.Watts = dt.conT0T1 %>% subset(Diet == "Watts")
# 
# pseq.ZIRC = ps.conT0T1 %>% subset_samples(Diet == "ZIRC")
# data.ZIRC = dt.conT0T1 %>% subset(Diet == "ZIRC")
# 
# 
# 
# # Genus level data
# Genus_data.Gemma = aggregate_taxa(pseq.Gemma, "Genus")
# # Family level data
# Family_data.Gemma = aggregate_taxa(pseq.Gemma, "Family")
# # Phylum level data
# Phylum_data.Gemma = aggregate_taxa(pseq.Gemma, "Phylum")
# 
# # Genus level data
# Genus_data.Watts = aggregate_taxa(pseq.Watts, "Genus")
# # Family level data
# Family_data.Watts = aggregate_taxa(pseq.Watts, "Family")
# # Phylum level data
# Phylum_data.Watts = aggregate_taxa(pseq.Watts, "Phylum")
# 
# # Genus level data
# Genus_data.ZIRC = aggregate_taxa(pseq.ZIRC, "Genus")
# # Family level data
# Family_data.ZIRC = aggregate_taxa(pseq.ZIRC, "Family")
# # Phylum level data
# Phylum_data.ZIRC = aggregate_taxa(pseq.ZIRC, "Phylum")
# 
# # Run Ancom
# 
# tax.label = "Genus"
# 
# tax.level.Gemma = Genus_data.Gemma 
# tax.data.Gemma = taxa.data.table(tax.level.Gemma)
# 
# tax.level.Watts = Genus_data.Watts 
# tax.data.Watts = taxa.data.table(tax.level.Watts)
# 
# tax.level.ZIRC = Genus_data.ZIRC 
# tax.data.ZIRC = taxa.data.table(tax.level.ZIRC)
# 
# # Genus_ Phylum_ Family_
# 
# # Gemma
# 
# out.Gemma = ancombc(phyloseq = tax.level.Gemma, 
#                     formula = "Timepoint", 
#                     p_adj_method = "BH", 
#                     # lib_cut = 0, 
#                     # prv_cut = 0,
#                     # group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#                     # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#                     neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#                     tol = 1e-05, 
#                     max_iter = 100,
#                     conserve = T, 
#                     alpha = 0.05, 
#                     # global = T
#                     #assay_name = "counts"
# )
# 
# res.Gemma = out.Gemma$res
# res_global.Gemma = out.Gemma$res_global
# 
# 
# # Watts
# 
# out.Watts = ancombc(phyloseq = tax.level.Watts, 
#                     formula = "Timepoint", 
#                     p_adj_method = "BH", 
#                     # lib_cut = 0, 
#                     # prv_cut = 0,
#                     # group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#                     # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#                     neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#                     tol = 1e-05, 
#                     max_iter = 100,
#                     conserve = T, 
#                     alpha = 0.05, 
#                     # global = T
#                     #assay_name = "counts"
# )
# 
# res.Watts = out.Watts$res
# res_global.Watts = out.Watts$res_global
# 
# 
# # ZIRC
# 
# out.ZIRC = ancombc(phyloseq = tax.level.ZIRC, 
#                    formula = "Timepoint", 
#                    p_adj_method = "BH", 
#                    # lib_cut = 0, 
#                    # prv_cut = 0,
#                    # group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#                    # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#                    neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#                    tol = 1e-05, 
#                    max_iter = 100,
#                    conserve = T, 
#                    alpha = 0.05, 
#                    # global = T
#                    #assay_name = "counts"
# )
# 
# res.ZIRC = out.ZIRC$res
# res_global.ZIRC = out.ZIRC$res_global
# 
# 
# 
# # Gemma
# 
# lfc.Gemma <- list()
# 
# lfc.Gemma$tab_lfc = res.Gemma$lfc
# lfc.Gemma$tab_se = res.Gemma$se
# lfc.Gemma$tab_p = res.Gemma$p_val
# lfc.Gemma$tab_diff = res.Gemma$diff_abn
# 
# lfc.Gemma$df_lfc = data.frame(lfc.Gemma$tab_lfc * lfc.Gemma$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# lfc.Gemma$df_se = data.frame(lfc.Gemma$tab_se * lfc.Gemma$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# colnames(lfc.Gemma$df_se)[-1] = paste0(colnames(lfc.Gemma$df_se)[-1], "SE")
# 
# lfc.Gemma$df_lfc <- 
#   lfc.Gemma$df_lfc %>%
#   rename(Gemma.6mpf = Timepoint6mpf)
# 
# lfc.Gemma$df_se <- 
#   lfc.Gemma$df_se %>%
#   rename(Gemma.6mpf = Timepoint6mpfSE)
# 
# # Watts
# 
# lfc.Watts <- list()
# 
# lfc.Watts$tab_lfc = res.Watts$lfc
# lfc.Watts$tab_se = res.Watts$se
# lfc.Watts$tab_p = res.Watts$p_val
# lfc.Watts$tab_diff = res.Watts$diff_abn
# 
# lfc.Watts$df_lfc = data.frame(lfc.Watts$tab_lfc * lfc.Watts$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# lfc.Watts$df_se = data.frame(lfc.Watts$tab_se * lfc.Watts$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# colnames(lfc.Watts$df_se)[-1] = paste0(colnames(lfc.Watts$df_se)[-1], "SE")
# 
# lfc.Watts$df_lfc <- 
#   lfc.Watts$df_lfc %>%
#   rename(Watts.6mpf = Timepoint6mpf)
# 
# lfc.Watts$df_se <- 
#   lfc.Watts$df_se %>%
#   rename(Watts.6mpf = Timepoint6mpfSE)
# 
# # ZIRC
# 
# lfc.ZIRC <- list()
# 
# lfc.ZIRC$tab_lfc = res.ZIRC$lfc
# lfc.ZIRC$tab_se = res.ZIRC$se
# lfc.ZIRC$tab_p = res.ZIRC$p_val
# lfc.ZIRC$tab_diff = res.ZIRC$diff_abn
# 
# lfc.ZIRC$df_lfc = data.frame(lfc.ZIRC$tab_lfc * lfc.ZIRC$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# lfc.ZIRC$df_se = data.frame(lfc.ZIRC$tab_se * lfc.ZIRC$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# colnames(lfc.ZIRC$df_se)[-1] = paste0(colnames(lfc.ZIRC$df_se)[-1], "SE")
# 
# lfc.ZIRC$df_lfc <- 
#   lfc.ZIRC$df_lfc %>%
#   rename(ZIRC.6mpf = Timepoint6mpf)
# 
# lfc.ZIRC$df_se <- 
#   lfc.ZIRC$df_se %>%
#   rename(ZIRC.6mpf = Timepoint6mpfSE)
# 
# # Statistical Data Tables
# lfc.time.diets <- list(lfc.Gemma,
#                        lfc.Watts,
#                        lfc.ZIRC) %>% 
#   setNames(c("Gemma", "Watts", "ZIRC"))
# 
# 
# # Join and Pivot
# 
# ## DF
# 
# df_lfc.dtime<- 
#   lfc.Gemma$df_lfc %>%
#   left_join(lfc.Watts$df_lfc, by = "Taxon_id",  na_matches = "never") %>% 
#   left_join(lfc.ZIRC$df_lfc, by = "Taxon_id",  na_matches = "never") %>% 
#   replace(is.na(.), 0) %>%
#   pivot_longer(cols = 2:4, names_to = "Diet")
# 
# ## SE
# 
# df_se.time <- 
#   lfc.Gemma$df_se %>%
#   left_join(lfc.Watts$df_se, by = "Taxon_id",  na_matches = "never") %>% 
#   left_join(lfc.ZIRC$df_se, by = "Taxon_id",  na_matches = "never") %>% 
#   replace(is.na(.), 0) %>%
#   pivot_longer(cols = 2:4, names_to = "Diet")
# 
# 
# # df_lfc.diets
# # df_se.diets
# 
# 
# 
# 
# 
# ## Diet:Exposure -----------------------------------------------------------
# 
# # Load Data
# 
# pseq.Gemma = ps.T1 %>% subset_samples(Diet == "Gemma")
# data.Gemma = dt.T1 %>% subset(Diet == "Gemma")
# 
# pseq.Watts = ps.T1 %>% subset_samples(Diet == "Watts")
# data.Watts = dt.T1 %>% subset(Diet == "Watts")
# 
# pseq.ZIRC = ps.T1 %>% subset_samples(Diet == "ZIRC")
# data.ZIRC = dt.T1 %>% subset(Diet == "ZIRC")
# 
# # Relevel Factors
# 
# levels(sample_data(pseq.Gemma)$PrePostExp) <- factor(levels(sample_data(pseq.Gemma)$PrePostExp),levels = c("Unexposed", "Exposed"))
# sample_data(pseq.Gemma)$PrePostExp <- relevel(factor(sample_data(pseq.Gemma)$PrePostExp), ref = "Unexposed")
# 
# levels(data.Gemma$PrePostExp) <- factor(levels(data.Gemma$PrePostExp), levels = c("Unexposed", "Exposed"))
# data.Gemma$PrePostExp <- relevel(factor(data.Gemma$PrePostExp), ref = "Unexposed")
# 
# levels(sample_data(pseq.Watts)$PrePostExp) <- factor(levels(sample_data(pseq.Watts)$PrePostExp),levels = c("Unexposed", "Exposed"))
# sample_data(pseq.Watts)$PrePostExp <- relevel(factor(sample_data(pseq.Watts)$PrePostExp), ref = "Unexposed")
# 
# levels(data.Watts$PrePostExp) <- factor(levels(data.Watts$PrePostExp), levels = c("Unexposed", "Exposed"))
# data.Watts$PrePostExp <- relevel(factor(data.Watts$PrePostExp), ref = "Unexposed")
# 
# levels(sample_data(pseq.ZIRC)$PrePostExp) <- factor(levels(sample_data(pseq.ZIRC)$PrePostExp),levels = c("Unexposed", "Exposed"))
# sample_data(pseq.ZIRC)$PrePostExp <- relevel(factor(sample_data(pseq.ZIRC)$PrePostExp), ref = "Unexposed")
# 
# levels(data.ZIRC$PrePostExp) <- factor(levels(data.ZIRC$PrePostExp), levels = c("Unexposed", "Exposed"))
# data.ZIRC$PrePostExp <- relevel(factor(data.ZIRC$PrePostExp), ref = "Unexposed")
# 
# 
# 
# # Genus level data
# Genus_data.Gemma = aggregate_taxa(pseq.Gemma, "Genus")
# # Family level data
# Family_data.Gemma = aggregate_taxa(pseq.Gemma, "Family")
# # Phylum level data
# Phylum_data.Gemma = aggregate_taxa(pseq.Gemma, "Phylum")
# 
# # Genus level data
# Genus_data.Watts = aggregate_taxa(pseq.Watts, "Genus")
# # Family level data
# Family_data.Watts = aggregate_taxa(pseq.Watts, "Family")
# # Phylum level data
# Phylum_data.Watts = aggregate_taxa(pseq.Watts, "Phylum")
# 
# # Genus level data
# Genus_data.ZIRC = aggregate_taxa(pseq.ZIRC, "Genus")
# # Family level data
# Family_data.ZIRC = aggregate_taxa(pseq.ZIRC, "Family")
# # Phylum level data
# Phylum_data.ZIRC = aggregate_taxa(pseq.ZIRC, "Phylum")
# 
# # Run Ancom
# 
# tax.label = "Genus"
# 
# tax.level.Gemma = Genus_data.Gemma 
# tax.data.Gemma = taxa.data.table(tax.level.Gemma)
# 
# tax.level.Watts = Genus_data.Watts 
# tax.data.Watts = taxa.data.table(tax.level.Watts)
# 
# tax.level.ZIRC = Genus_data.ZIRC 
# tax.data.ZIRC = taxa.data.table(tax.level.ZIRC)
# 
# # Genus_ Phylum_ Family_
# 
# # Gemma
# 
# out.Gemma = ancombc(phyloseq = tax.level.Gemma, 
#                     formula = "PrePostExp", 
#                     p_adj_method = "BH", 
#                     # lib_cut = 0, 
#                     # prv_cut = 0,
#                     # group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#                     # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#                     neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#                     tol = 1e-05, 
#                     max_iter = 100,
#                     conserve = T, 
#                     alpha = 0.05, 
#                     # global = T
#                     #assay_name = "counts"
# )
# 
# res.Gemma = out.Gemma$res
# res_global.Gemma = out.Gemma$res_global
# 
# 
# # Watts
# 
# out.Watts = ancombc(phyloseq = tax.level.Watts, 
#                     formula = "PrePostExp", 
#                     p_adj_method = "BH", 
#                     # lib_cut = 0, 
#                     # prv_cut = 0,
#                     # group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#                     # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#                     neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#                     tol = 1e-05, 
#                     max_iter = 100,
#                     conserve = T, 
#                     alpha = 0.05, 
#                     # global = T
#                     #assay_name = "counts"
# )
# 
# res.Watts = out.Watts$res
# res_global.Watts = out.Watts$res_global
# 
# 
# # ZIRC
# 
# out.ZIRC = ancombc(phyloseq = tax.level.ZIRC, 
#                    formula = "PrePostExp", 
#                    p_adj_method = "BH", 
#                    # lib_cut = 0, 
#                    # prv_cut = 0,
#                    # group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
#                    # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
#                    neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
#                    tol = 1e-05, 
#                    max_iter = 100,
#                    conserve = T, 
#                    alpha = 0.05, 
#                    # global = T
#                    #assay_name = "counts"
# )
# 
# res.ZIRC = out.ZIRC$res
# res_global.ZIRC = out.ZIRC$res_global
# 
# # Gemma
# 
# lfc.Gemma <- list()
# 
# lfc.Gemma$tab_lfc = res.Gemma$lfc
# lfc.Gemma$tab_se = res.Gemma$se
# lfc.Gemma$tab_p = res.Gemma$p_val
# lfc.Gemma$tab_diff = res.Gemma$diff_abn
# 
# lfc.Gemma$df_lfc = data.frame(lfc.Gemma$tab_lfc * lfc.Gemma$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# lfc.Gemma$df_se = data.frame(lfc.Gemma$tab_se * lfc.Gemma$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# colnames(lfc.Gemma$df_se)[-1] = paste0(colnames(lfc.Gemma$df_se)[-1], "SE")
# 
# lfc.Gemma$df_lfc <- 
#   lfc.Gemma$df_lfc %>%
#   rename(Gemma.Exposed = PrePostExpExposed)
# 
# lfc.Gemma$df_se <- 
#   lfc.Gemma$df_se %>%
#   rename(Gemma.Exposed = PrePostExpExposedSE)
# 
# # Watts
# 
# lfc.Watts <- list()
# 
# lfc.Watts$tab_lfc = res.Watts$lfc
# lfc.Watts$tab_se = res.Watts$se
# lfc.Watts$tab_p = res.Watts$p_val
# lfc.Watts$tab_diff = res.Watts$diff_abn
# 
# lfc.Watts$df_lfc = data.frame(lfc.Watts$tab_lfc * lfc.Watts$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# lfc.Watts$df_se = data.frame(lfc.Watts$tab_se * lfc.Watts$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# colnames(lfc.Watts$df_se)[-1] = paste0(colnames(lfc.Watts$df_se)[-1], "SE")
# 
# lfc.Watts$df_lfc <- 
#   lfc.Watts$df_lfc %>%
#   rename(Watts.Exposed = PrePostExpExposed)
# 
# lfc.Watts$df_se <- 
#   lfc.Watts$df_se %>%
#   rename(Watts.Exposed = PrePostExpExposedSE)
# 
# # ZIRC
# 
# lfc.ZIRC <- list()
# 
# lfc.ZIRC$tab_lfc = res.ZIRC$lfc
# lfc.ZIRC$tab_se = res.ZIRC$se
# lfc.ZIRC$tab_p = res.ZIRC$p_val
# lfc.ZIRC$tab_diff = res.ZIRC$diff_abn
# 
# lfc.ZIRC$df_lfc = data.frame(lfc.ZIRC$tab_lfc * lfc.ZIRC$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# lfc.ZIRC$df_se = data.frame(lfc.ZIRC$tab_se * lfc.ZIRC$tab_diff, check.names = FALSE) %>% 
#   rownames_to_column("Taxon_id")
# colnames(lfc.ZIRC$df_se)[-1] = paste0(colnames(lfc.ZIRC$df_se)[-1], "SE")
# 
# lfc.ZIRC$df_lfc <- 
#   lfc.ZIRC$df_lfc %>%
#   rename(ZIRC.Exposed = PrePostExpExposed)
# 
# lfc.ZIRC$df_se <- 
#   lfc.ZIRC$df_se %>%
#   rename(ZIRC.Exposed = PrePostExpExposedSE)
# 
# # Statistical Data Tables
# lfc.exp.diets <- list(lfc.Gemma,
#                       lfc.Watts,
#                       lfc.ZIRC) %>% 
#   setNames(c("Gemma", "Watts", "ZIRC"))
# 
# 
# # Join and Pivot
# 
# ## DF
# 
# df_lfc.exp<- 
#   lfc.Gemma$df_lfc %>%
#   left_join(lfc.Watts$df_lfc, by = "Taxon_id",  na_matches = "never") %>% 
#   left_join(lfc.ZIRC$df_lfc, by = "Taxon_id",  na_matches = "never") %>% 
#   replace(is.na(.), 0) %>%
#   pivot_longer(cols = 2:4, names_to = "Diet")
# 
# ## SE
# 
# df_se.exp <- 
#   lfc.Gemma$df_se %>%
#   left_join(lfc.Watts$df_se, by = "Taxon_id",  na_matches = "never") %>% 
#   left_join(lfc.ZIRC$df_se, by = "Taxon_id",  na_matches = "never") %>% 
#   replace(is.na(.), 0) %>%
#   pivot_longer(cols = 2:4, names_to = "Diet")
# 
# 
# 
# 
# 
# 
# 

# End of doc ---------------------------------------------------------------------

save.image(file.path(path.data, 
                     paste0("/R_objects/RoL_HeaterTrial_PostMicrobiomeProcessing_ENV_", 
                            Sys.Date(),
                            ".RData"
                            ) ) )
