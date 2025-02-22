###binary_models.R
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

data$Name <- as.factor(data$Name)
data$Sex <- as.factor(data$Sex)
data$AgeClass <- as.factor(data$AgeClass)
data$ActiveFungus <- as.factor(data$ActiveFungus)
data$Name2 <- data$Name
data$Name3 <- data$Name

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

#Set Priors
priorl1 <- set_prior("normal(0,10)", class = "b") +
  set_prior("cauchy(0,5)", class = "sd") #+ 
# set_prior("gamma(0.01,0.01)", class ="shape")

#Formula
f.ID_gm_hro <- bf(ActiveFungus ~ Sex + Time.in.pop + InfectedDensity.Annual + 
            (1 | gr(Name, cov = gm2_hr)) + (1 | gr(Name2, cov = pd_hrm)) + (1 | Name3) + (1 | Dragon.Year))
f.gm_hro <- bf(ActiveFungus ~ Sex + Time.in.pop + 
            (1 | gr(Name, cov = gm2_hr)) + (1 | gr(Name2, cov = pd_hrm)) + (1 | Name3) + (1 | Dragon.Year))
f.ID_gm <- bf(ActiveFungus ~ Sex + Time.in.pop + InfectedDensity.Annual + 
                (1 | gr(Name, cov = gm2)) + (1 | Name2) + (1 | Dragon.Year))
f.gm <- bf(ActiveFungus ~ Sex + Time.in.pop + 
             (1 | gr(Name, cov = gm2)) + (1 | Name2) + (1 | Dragon.Year))

f.1l <- bf(InfectedDensity.Annual.log ~ 1 + Sex + Time.in.pop + ActiveFungus + 
             (1 | gr(Name, cov = gm2)) + (1 | Name2) + (1 | Dragon.Year))
f.1g_hro <- bf(InfectedDensity.Annual.adj ~ 1 + Sex + Time.in.pop + ActiveFungus + (1 | gr(Name, cov = pd_hrm)) + (1 | Name2) + (1 | Dragon.Year))


#Run Model
AF_ID_grm_hro <- brm(
  f.ID_gm_hro,
  prior = priorl1,
  data = data_HRO,
  data2 = list(gm2_hr = gm2_hr, pd_hrm = pd_hrm),
  family = bernoulli(),
  chains = 4, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Binary/brm_ID_GRM_HRO_V2"
)

summary(AF_ID_grm_hro)
plot(AF_ID_grm_hro)

AF_noID_grm_hro <- brm(
  f.gm_hro,
  prior = priorl1,
  data = data_HRO,
  data2 = list(gm2_hr = gm2_hr, pd_hrm = pd_hrm),
  family = bernoulli(),
  chains = 4, cores = 2, iter = 15000, thin = 4, warmup = 5000, 
  threads = threading(2),
  set.seed(2595), file = "Output/Binary/brm_GRM_HRO_V2"
)

AF_ID_grm <- brm(
  f.ID_gm,
  prior = priorl1,
  data = data_GRM,
  data2 = list(gm2 = gm2),
  family = bernoulli(),
  chains = 4, cores = 2, iter = 15000, thin = 4, warmup = 5000,
  threads = threading(2),
  set.seed(2595), file = "Output/Binary/brm_ID_GRM_V1"
)

summary(AF_ID_grm_hro)
plot(AF_ID_grm_hro)

AF_noID_grm <- brm(
  f.gm,
  prior = priorl1,
  data = data_GRM,
  data2 = list(gm2 = gm2),
  family = bernoulli(),
  chains = 4, cores = 2, iter = 15000, thin = 4, warmup = 5000, 
  threads = threading(2),
  set.seed(2595), file = "Output/Binary/brm_GRM_V1"
)

summary(AF_noID_grm_hro)
plot(AF_noID_grm_hro)

loo(AF_noID_grm_hro)
loo(AF_ID_grm_hro)

n_parameters(AF_ID_grm_hro)
mcmc_plot(AF_ID_grm_hro, type = "acf")

AF_ID_grm_hro <- add_criterion(AF_ID_grm_hro, "loo")
AF_noID_grm_hro <- add_criterion(AF_noID_grm_hro, "loo")
AF_ID_grm <- add_criterion(AF_ID_grm, "loo")
AF_noID_grm <- add_criterion(AF_noID_grm, "loo")

loo_compare(AF_ID_grm_hro, AF_noID_grm_hro, criterion = "loo")
loo_compare(AF_ID_grm, AF_noID_grm, criterion = "loo")

Var.table <- as_draws_df(AF_ID_grm_hro)

Var.table$h.bwt.1 <- as.mcmc((Var.table$sd_Name__Intercept)^2 / 
                               ((Var.table$sd_Name__Intercept)^2 + 
                                  (Var.table$sd_Name2__Intercept)^2 + 
                                  (Var.table$sd_Dragon.Year__Intercept)^2 + 1))

summary(Var.table$h.bwt.1)
plot(Var.table$h.bwt.1)

yhat <- predict(AF_ID, type = "terms")

pr <-  purrr::pmap_dfr(list(mu = Var.table$Intercept, # because i have fixed effects do i not need to average over them?
                            var.a = Var.table$sd_Name__Intercept, 
                            var.p = rowSums(Var.table[5:7]) + 1), #need to add the residual 
                       QGparams, custom.model = "binom1.probit", verbose = F) 
mean(pr[["h2.obs"]]) 
HPDinterval(as.mcmc(pr[["h2.obs"]]))



