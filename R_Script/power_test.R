###power_analysis.R
library(brms); library(parallel); library(corpcor); library(lubridate); library(tidyverse)
library(magrittr); library(QGglmm); library(purrr); library(MCMCglmm); library(insight); library(tidybayes)

if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
    Sys.info()["sysname"] == "MacOS" && getRversion() == "4.4.2") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

# #Making GRM
# g2<-read.csv("Data/R_estimates_drag.csv" ,stringsAsFactors = FALSE)
# table1<-aggregate(g2$ind1.id,list(g2$ind1.id),length)
# table2<-aggregate(g2$ind2.id,list(g2$ind2.id),length)
# names(table1)<-c("Name", "Number")
# table1<-table1[order(table1$Number,decreasing=TRUE),]
# names(table2)<-c("Name", "Number")
# table2<-table2[order(table2$Number,decreasing=FALSE),]
# g2b<-g2[order(match(g2$ind2.id, table2$Name)),]
# g2b<-g2b[order(match(g2b$ind1.id, table1$Name)),]
# nam=unique(c(as.character(table1$Name),as.character(table2$Name)))
# gm<-matrix(0,880,880)
# rownames(gm)<-colnames(gm)<-nam
# gm[lower.tri(gm,diag=FALSE)]<-g2b$quellergt
# diag(gm)<-1.001
# gm[upper.tri(gm)]<-t(gm)[upper.tri(gm)]
# gm <- gm[!(rownames(gm) %in% name_to_rm),
#          !(colnames(gm) %in% name_to_rm)]
# indsm <- rownames(gm)

#GRM
load("Data/pd_GRM.RData")

#HROM
load("Data/pd_HRO95_3.RData")

#Pheno Data
data <- read.csv("Data/One_per_season.csv", header = T)

names_to_keep <- intersect(rownames(pd_hrm), data$Name)
data <- data[data$Name %in% names_to_keep, ]
pd_hrm <- pd_hrm[rownames(pd_hrm) %in% data_HRO$Name,
                 colnames(pd_hrm) %in% data_HRO$Name]
gm2_hr <- gm2[rownames(gm2) %in% data_HRO$Name,
              colnames(gm2) %in% data_HRO$Name]
correlation <- cor(c(gm2_hr), c(pd_hrm), method = "pearson")

data$genetic <- rnorm(nrow(data), mean = 0, sd = 1) 
data$genetic <- gm2_hr %*% data$genetic 

data$environment <- rnorm(nrow(data), mean = 0, sd = 0.5) 
data$environment <- pd_hrm %*% data$environment 

phenotype <- genetic + envrironment 

#Combine Similarity Matrices
weight_genetic <- 0.7  
weight_space <- 0.3
combined_sim_matrix <- (weight_genetic * gm2_hr + weight_space * pd_hrm)

# Convert similarity matrices to dissimilarity matrices
dissimilarity_genetic <- 1 - gm2_hr
dissimilarity_space <- 1 - pd_hrm

# Perform Mantel test to assess the relationship between the two dissimilarity matrices
mantel_result <- mantel(dissimilarity_genetic, dissimilarity_space, permutations = 9999) # Increased permutations for more robust p-value

# Extract residuals from the Mantel test (these represent the variation not explained by the relationship between the two)
mantel_residuals <- mantel_result$residual

# Convert the residuals back to a similarity matrix. The logic here is that
# smaller absolute residuals indicate a stronger combined similarity. We use a
# Gaussian kernel to convert the residuals into similarities. You might need to
# adjust the bandwidth parameter (sigma) depending on the scale of your residuals.
# A good starting point is the standard deviation of your residuals.

sigma <- sd(mantel_residuals) # Bandwidth parameter (adjust as needed)
combined_sim_matrix <- matrix(0, nrow = length(common_individuals), ncol = length(common_individuals), dimnames = list(common_individuals, common_individuals))

for (i in 1:length(common_individuals)) {
  for (j in 1:length(common_individuals)) {
    combined_sim_matrix[i, j] <- exp(-((mantel_residuals[i, j])^2) / (2 * sigma^2))
  }
}

# 3. Create the New Variable in the Dataframe

# Initialize the new column
data$combined_similarity_score <- NA
i=1
for (i in 1:nrow(data)) {
  current_individual <- data$Name[i]
  
  # Calculate the similarity score for the current individual
  # by comparing it to all other individuals *within the current year*
  current_year <- data$Dragon.Year[i]
  individuals_in_year <- data$Name[ data$Dragon.Year == current_year]
  
  similarities <- combined_sim_matrix[current_individual, individuals_in_year]
  
  # Choose how to summarize similarity (mean, median, max, etc.)
  data$combined_similarity_score[i] <- mean(similarities, na.rm=TRUE)
  
  #Or if you want to take into account all other individuals across all years:
  #similarities <- combined_sim_matrix[current_individual, ]
  #df$combined_similarity_score[i] <- mean(similarities, na.rm=TRUE)
}

# 4. (Optional) Validate using Mantel test (not necessary for this approach, but can be informative)
# You can compare the combined similarity to the original ones to see how well they correlate.

dissimilarity_combined <- 1 - combined_sim_matrix # Convert to dissimilarities
mantel(dissimilarity_combined, dissimilarity_genetic)
mantel(dissimilarity_combined, dissimilarity_space)


# Print or inspect the updated dataframe
print(head(df))

data %>% ggplot(aes(x = sim_trait)) + geom_density()

data$Name <- as.factor(data$Name)
data$Sex <- as.factor(data$Sex)
data$AgeClass <- as.factor(data$AgeClass)
data$ActiveFungus <- as.factor(data$ActiveFungus)
data$Name2 <- data$Name
data$Name3 <- data$Name

Classes <- data %>% dplyr::select(Time.in.pop, InfectedDensity.Annual, sim_trait) %>% 
  sapply(class)
ToScale <- names(Classes[Classes %in% c("integer", "numeric")])
data[, paste0(ToScale, ".Original")] <- data[, ToScale]
data %<>% mutate_at(ToScale, ~c(scale(.x)))

data <- data %>%
  group_by(Dragon.Year) %>%
  mutate(InfectedDensity.Relative = InfectedDensity.Annual.Original - mean(InfectedDensity.Annual.Original),
         InfectedDensity.Normalized = scale(InfectedDensity.Annual.Original)) %>%
  ungroup()

#Set Priors
priorl1 <- set_prior("normal(0,10)", class = "b") +
  set_prior("cauchy(0,5)", class = "sd") #+ 
# set_prior("gamma(0.01,0.01)", class ="shape")

f.pt_grm_hro <- bf(combined_similarity_score ~ Sex + Time.in.pop + ActiveFungus + 
            (1 | gr(Name, cov = gm2_hr)) + (1 | gr(Name2, cov = pd_hrm)) + (1 | Name3) + (1 | Dragon.Year))

f.pt_hro <- bf(combined_similarity_score ~ Sex + Time.in.pop + ActiveFungus + 
                     (1 | gr(Name, cov = pd_hrm)) + (1 | Name2) + (1 | Dragon.Year))

f.pt_grm <- bf(combined_similarity_score ~ Sex + Time.in.pop + ActiveFungus + 
                     (1 | gr(Name, cov = gm2_hr)) + (1 | Name2) + (1 | Dragon.Year))

power_test_grm_hro <- brm(
  f.pt_grm_hro,
  prior = priorl1,
  data = data,
  data2 = list(gm2_hr = gm2_hr, pd_hrm = pd_hrm),
  # family = Gamma(link = "log"),
  chains = 2, cores = 2, iter = 30000, thin = 4, warmup = 9000,
  threads = threading(2),
  set.seed(2595), file = "Output/brms_powertest_V4"
)

power_test_hro <- brm(
  f.pt_hro,
  prior = priorl1,
  data = data,
  data2 = list(pd_hrm = pd_hrm),
  # family = Gamma(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/brms_powertest_hro_V1"
)


power_test_grm <- brm(
  f.pt_grm,
  prior = priorl1,
  data = data,
  data2 = list(gm2_hr = gm2_hr),
  # family = Gamma(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/brms_powertest_grm_V1"
)

plot(power_test_grm_hro)
summary(power_test_grm_hro)
loo(power_test)

res <- get_variance_residual(power_test_grm)

Var.table <- as_draws_df(power_test_grm)

Var.table$h.bwt.1 <- as.mcmc((Var.table$sd_Name__Intercept)^2 / 
                               ((Var.table$sd_Name__Intercept)^2 + 
                                  (Var.table$sd_Name2__Intercept)^2 + 
                                  (Var.table$sd_Dragon.Year__Intercept)^2))

summary(Var.table$h.bwt.1)
plot(Var.table$h.bwt.1)


# Function to assign weighted random values
assign_weighted_random <- function(data, similarity_matrix) {
  n_individuals <- nrow(data)
  weighted_random_values <- numeric(n_individuals)
  
  for (i in 1:n_individuals) {
    current_name <- data$Name[i]
    similar_names <- names(gm2_hr[current_name,][order(gm2_hr[current_name, ], decreasing = TRUE)]) #Order by similarity
    similar_names <- similar_names[similar_names != current_name] #Remove the current name from the similar names
    
    weights <- gm2_hr[current_name, similar_names] #Extract the similarity weights
    weights <- weights/sum(weights) #Normalise the weights so they sum to 1
    
    random_values <- runif(length(similar_names)) #Generate random numbers for the similar names
    weighted_values <- random_values * weights #Weight the random numbers by similarity
    
    #Option 1: Average the weighted values of similar individuals
    weighted_random_values[i] <- mean(weighted_values)
    
    #Option 2: Sample from similar individuals' values (if they already exist)
    #if ("weighted_random" %in% colnames(df)) { #Check if the column exists
    #  existing_values <- df$weighted_random[df$Name %in% similar_names]
    #  if (length(existing_values) > 0) { #Check if there are values to sample from
    #    weighted_random_values[i] <- sample(existing_values, 1, prob = weights)
    #  } else {
    #    weighted_random_values[i] <- mean(weighted_values) #Fallback to averaging
    #  }
    #} else {
    #  weighted_random_values[i] <- mean(weighted_values) #Fallback to averaging
    #}
  }
  
  data$sim_trait <- weighted_random_values
  return(data)
}


df <- assign_weighted_random(df, similarity_matrix)
print(df)

