###Years_diseased_healthy.R
library(brms); library(parallel); library(corpcor); library(lubridate); library(tidyverse)
library(magrittr); library(QGglmm); library(purrr); library(MCMCglmm); library(insight): library(tidybayes)

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

data$Name <- as.factor(data$Name)
data$Sex <- as.factor(data$Sex)
data$AgeClass <- as.factor(data$AgeClass)
data$ActiveFungus <- as.factor(data$ActiveFungus)
data$Name2 <- data$Name
data$Name3 <- data$Name

data <- data %>%
  arrange(Name, Dragon.Year) %>%
  group_by(Name, Dragon.Year) %>%  
  mutate(
    ActiveFungus_Final = max(ActiveFungus) 
  ) %>%
  ungroup() %>%
  group_by(Name) %>%
  mutate(
    Years.healthy = {
      years_healthy <- rep(0, n())
      cumulative_healthy_years <- 0
      
      for (i in 1:n()) {
        if (ActiveFungus_Final[i] == 0) {
          cumulative_healthy_years <- cumulative_healthy_years + 1
          years_healthy[i] <- cumulative_healthy_years
        } else {
          cumulative_healthy_years <- 0
          years_healthy[i] <- 0
        }
      }
      years_healthy
    }
  ) %>%
  ungroup() %>%
  select(-ActiveFungus_Final) 

data <- data %>%
  arrange(Name, Dragon.Year) %>%
  group_by(Name, Dragon.Year) %>%  
  mutate(
    ActiveFungus_Final = max(ActiveFungus) 
  ) %>%
  ungroup() %>%
  group_by(Name) %>% 
  mutate(
    Years.diseased = {
      years_diseased <- rep(0, n())
      cumulative_diseased_years <- 0
      
      for (i in 1:n()) {
        if (ActiveFungus_Final[i] == 1) {
          cumulative_diseased_years <- cumulative_diseased_years + 1
          years_diseased[i] <- cumulative_diseased_years
        } else {
          cumulative_diseased_years <- 0
          years_diseased[i] <- 0
        }
      }
      years_diseased
    }
  ) %>%
  ungroup() %>%
  select(-ActiveFungus_Final)

data <- data %>%
  arrange(Name, Dragon.Year) %>%
  group_by(Name, Dragon.Year) %>%
  mutate(
    ActiveFungus_Final = max(ActiveFungus)
  ) %>%
  group_by(Name) %>%
  mutate(
    Years.diseased = cumsum(ActiveFungus_Final),
    Years.diseased = ifelse(ActiveFungus_Final == 0, 0, Years.diseased),  # Reset to 0 for healthy years
    Years.healthy = {  # Calculate Years.healthy
      years_healthy <- rep(0, n())
      cumulative_healthy_years <- 0
      
      for (i in 1:n()) {
        if (ActiveFungus_Final[i] == 0) {
          cumulative_healthy_years <- cumulative_healthy_years + 1
          years_healthy[i] <- cumulative_healthy_years
        } else {
          cumulative_healthy_years <- 0
          years_healthy[i] <- 0
        }
      }
      years_healthy
    }
  ) %>%
  ungroup() #%>%
  select(-ActiveFungus_Final)

data %>% ggplot(aes(x = Years.diseased)) + geom_bar() 
data %>% ggplot(aes(x = Years.healthy)) + geom_bar() #+ facet_wrap(data$Dragon.Year)

Classes <- data %>% dplyr::select(Time.in.pop, InfectedDensity.Annual) %>% 
  sapply(class)
ToScale <- names(Classes[Classes %in% c("integer", "numeric")])
data[, paste0(ToScale, ".Original")] <- data[, ToScale]
data %<>% mutate_at(ToScale, ~c(scale(.x)))

names_to_keep <- intersect(rownames(pd_hrm), data$Name)
data_HRO <- data[data$Name %in% names_to_keep, ]
data_GRM <- data[data$Name %in% rownames(gm2), ]
pd_hrm <- pd_hrm[rownames(pd_hrm) %in% data_HRO$Name,
                 colnames(pd_hrm) %in% data_HRO$Name]
gm2_hr <- gm2[rownames(gm2) %in% data_HRO$Name,
              colnames(gm2) %in% data_HRO$Name]
gm2 <- gm2[rownames(gm2) %in% data_GRM$Name,
           colnames(gm2) %in% data_GRM$Name]
correlation <- cor(c(gm2_hr), c(pd_hrm), method = "pearson")

priorl1 <- set_prior("normal(0,10)", class = "b") +
  set_prior("cauchy(0,5)", class = "sd") #+ 
# set_prior("gamma(0.01,0.01)", class ="shape")

f.YH_ID_grm_hro <- bf(Years.healthy ~ Sex + Time.in.pop + InfectedDensity.Annual + 
            (1 | gr(Name, cov = gm2_hr)) + (1 | gr(Name2, cov = pd_hrm)) + (1 | Name3) + (1 | Dragon.Year))

f.YH_ID_grm <- bf(Years.healthy ~ Sex + Time.in.pop + InfectedDensity.Annual + 
             (1 | gr(Name, cov = gm2))  + (1 | Name2) + (1 | Dragon.Year))

f.YH_grm <- bf(Years.healthy ~ Sex + Time.in.pop + 
                 (1 | gr(Name, cov = gm2))  + (1 | Name2) + (1 | Dragon.Year))

f.YD_ID_grm_hro <- bf(Years.diseased ~ Sex + Time.in.pop + InfectedDensity.Annual + 
                     (1 | gr(Name, cov = gm2_hr)) + (1 | gr(Name2, cov = pd_hrm)) + (1 | Name3) + (1 | Dragon.Year))

f.YD_ID_grm <- bf(Years.diseased ~ Sex + Time.in.pop + InfectedDensity.Annual + 
                 (1 | gr(Name, cov = gm2))  + (1 | Name2) + (1 | Dragon.Year))

f.YD_grm <- bf(Years.diseased ~ Sex + Time.in.pop + 
                 (1 | gr(Name, cov = gm2))  + (1 | Name2) + (1 | Dragon.Year))

YearsHealthy_ID_grm_hro <- brm(
  f.YH_ID_grm_hro,
  prior = priorl1,
  data = data_HRO,
  data2 = list(gm2_hr = gm2_hr, pd_hrm = pd_hrm),
  family = poisson(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Poisson/brm_YH_ID_GRM_HRO_V1"
)

YearsHealthy_ID_grm <- brm(
  f.YH_ID_grm,
  prior = priorl1,
  data = data_GRM,
  data2 = list(gm2 = gm2),
  family = poisson(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Poisson/brm_YH_ID_GRM_V1"
)

YearsHealthy_grm <- brm(
  f.YH_grm,
  prior = priorl1,
  data = data_GRM,
  data2 = list(gm2 = gm2),
  family = poisson(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Poisson/brm_YH_GRM_V2"
)

YearsDiseased_ID_grm_hro <- brm(
  f.YD_ID_grm_hro,
  prior = priorl1,
  data = data_HRO,
  data2 = list(gm2_hr = gm2_hr, pd_hrm = pd_hrm),
  family = poisson(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Poisson/brm_YD_ID_GRM_HRO_V1"
)

YearsDiseased_ID_grm <- brm(
  f.YD_ID_grm,
  prior = priorl1,
  data = data_GRM,
  data2 = list(gm2 = gm2),
  family = poisson(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Poisson/brm_YD_ID_GRM_V1"
)

YearsDiseased_grm <- brm(
  f.YD_grm,
  prior = priorl1,
  data = data_GRM,
  data2 = list(gm2 = gm2),
  family = poisson(link = "log"),
  chains = 2, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Poisson/brm_YD_GRM_V1"
)

summary(YearsHealthy_grm)
plot(YearsHealthy_ID_grm)

loo(YearsHealthy2)
