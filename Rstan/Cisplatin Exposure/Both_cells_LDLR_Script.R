library(conflicted)
conflict_prefer("filter","dplyr")
conflict_prefer("lag","dplyr")
library(cmdstanr)
library(bayesplot)
library(tidyverse)
library(deSolve)
library(rstudioapi)

options(mc.cores = parallel::detectCores())

setwd(dirname(getActiveDocumentContext()$path))



SE_calc <- function(x){(sd(x)/sqrt(n()))/5}


Score_data_CPT = read.csv("Inputdata/Cisplatin_Exposure_data.csv") %>% filter(Segment == "CPT") %>% group_by(timepoints,StateVar) %>% summarize(Score=mean(data4modelInterpol)) %>% filter(timepoints !=1)

SD_data_CPT = read.csv("Inputdata/Cisplatin_Exposure_data.csv") %>% filter(Segment == "CPT") %>% group_by(timepoints,StateVar) %>% summarize(SD= SE_calc(data4modelInterpol)) %>% filter(timepoints !=1)

DDRep_data_CPT <- Score_data_CPT %>% filter(StateVar == "DD" | StateVar == "Rep") %>% pivot_wider(names_from = "StateVar", values_from = "Score") %>% ungroup()#  %>% select(-timepoints) %>% as.matrix()

DDRep_sd_CPT <- SD_data_CPT %>% filter(StateVar == "DD" | StateVar == "Rep") %>% pivot_wider(names_from = "StateVar", values_from = "SD") %>% ungroup()#  %>% select(-timepoints) %>% as.matrix()

timepoints_data_CPT <- Score_data_CPT %>% filter(StateVar == "DD" | StateVar == "Rep") %>% pivot_wider(names_from = "StateVar", values_from = "Score") %>% ungroup() %>% pull(timepoints)


Score_data_PPT = read.csv("Inputdata/Cisplatin_Exposure_data.csv") %>% filter(Segment == "PPT") %>% group_by(timepoints,StateVar) %>% arrange(.by_group = TRUE) %>% summarize(Score=mean(data4modelInterpol)) %>% filter(timepoints !=1& timepoints !=120 & timepoints !=192)

SD_data_PPT = read.csv("Inputdata/Cisplatin_Exposure_data.csv") %>% filter(Segment == "PPT") %>% group_by(timepoints,StateVar) %>% summarize(SD= SE_calc(data4modelInterpol)) %>% filter(timepoints !=1 &timepoints !=120 & timepoints !=192)

DDRep_data_PPT <- Score_data_PPT %>% filter(StateVar == "DD" | StateVar == "Rep") %>% pivot_wider(names_from = "StateVar", values_from = "Score") %>% ungroup()#  %>% select(-timepoints) %>% as.matrix()

DDRep_sd_PPT <- SD_data_PPT %>% filter(StateVar == "DD" | StateVar == "Rep") %>% pivot_wider(names_from = "StateVar", values_from = "SD") %>% ungroup() #  %>% select(-timepoints) %>% as.matrix()

#DDRep_data_PPT[3,2] <- 1.25

timepoints_data_PPT <- Score_data_PPT %>% filter(StateVar == "DD" | StateVar == "Rep") %>% pivot_wider(names_from = "StateVar", values_from = "Score") %>% ungroup() %>% pull(timepoints)



data_list = list(
  N_CPT = length(timepoints_data_CPT),
  t0_CPT = 0,
  ts_CPT = unique(timepoints_data_CPT),
  y1= DDRep_data_CPT%>% select(-timepoints) %>% as.matrix(),
  y0_CPT= c(50,0,0,0,0),
  sigma1 = DDRep_sd_CPT%>% select(-timepoints) %>% as.matrix(),
  
  N_PPT = length(timepoints_data_PPT),
  t0_PPT = 0,
  ts_PPT = unique(timepoints_data_PPT),
  y2= DDRep_data_PPT%>% select(-timepoints) %>% as.matrix(),
  y0_PPT= c(50,0,0,0,0),
  sigma2 = DDRep_sd_PPT%>% select(-timepoints) %>% as.matrix())



options(mc.cores = parallel::detectCores())  
chains <- 4
iter <- 10000
warmup <- 5000


name_model <- "Both_Cells_LDLR"
compiled_model <- cmdstan_model(paste(name_model,".stan",sep=""))

Both_Cells_LDMR<- compiled_model$sample(data = data_list, 
                                         chains = chains,
                                         parallel_chains = getOption("mc.cores", 1), 
                                         iter_warmup = warmup, 
                                         iter_sampling = iter-warmup,
                                         show_messages = TRUE,
                                         refresh = 20,
                                         save_warmup = TRUE)

Both_Cells_LDMR$save_object(file = paste("Output/Full_",iter,name_model,format(Sys.Date()),".RDS",sep=""))
