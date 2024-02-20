#//////////////////////////////////////////////////////////////////////////////////////
# Aim: Create estimations of TOTs with SE using bootstapping
# Date: 20-februarr-2024
# Author: Javier Amaya-Nieto 
# Project: Cost-Sharing in Medical Care Can Increase Adult Mortality: Evidence from Colombia
#//////////////////////////////////////////////////////////////////////////////////////



# Montando las librerias
library(tidyverse)
library(ggplot2)
library(haven)
library(dplyr)
library(readxl)
library(magrittr)
library(reshape2)
library(rio)
library(openxlsx)
library(haven)
library(gridExtra)
library(extrafont)
font_import()
loadfonts(device = "win")


#### SETTING PARAMETERS CODE ####
results_path <- "D:/Documents_2/......."
files_path <- "D:/Documents_2/........."

## Establishing important variables for running the code
# Number of resamples performed
num_sam <- 500
num_lags <- 96 #INCLUDING LAG 0


## Defining the name of variables for y labels
name_y_label_list <- list("Outpatient Services")

## Creating a list of the variables to include
mysheetlist <- list("AMB")

#### SE using bootstrap ####

## indicator to move the label in the graphs
ind_ylabel_graphs <- 1 

for(i in mysheetlist) {
# The code run the estimation for each replication and then collect a matrix of the TOT calculated to lately, get the SE


  # Setting matrices to collect results

  mat_tot = matrix(ncol = num_sam, nrow = num_lags)
  BD_TOT_is <- data.frame(mat_tot)
  names(BD_TOT_is) <- paste("TOT_i", 1:num_sam, sep="")
  
  # Cumulative TOTs
  BD_TOT_cs <- data.frame(mat_tot)
  names(BD_TOT_cs) <- paste("TOT_c", 1:num_sam, sep="")
  
  # Final TOTs with SE
  mat_tot_f = matrix(ncol = 9, nrow = num_lags)
  BD_tot_se <- data.frame(mat_tot_f)
  names(BD_tot_se) <- c("Lag","TOTc","tot_SE","CI_low", "CI_upp", "ITT","itt_SE","CI_low_itt","CI_upp_itt")
  
  ind_itt <- 1
  
  while(ind_itt < num_sam+1){
    
    # Loaing data from TOTs and PIs
    df_tot <- read_dta(paste0(files_path,
                              "/Bootstrap/PI_",i,"_matrix_bootstrapping_LLR_500.dta"))
    
    
    ###////////////////////////////////////////////////////////////////////////
    # Making the code to repeat the estimations for each iteration of the bootstrap
    Lag_data <- df_tot[,1]
    ITT_data <- df_tot[,ind_itt+1] # indicador movil #1
    PI_data <- df_tot[,ind_itt+(num_sam+1)] # indicador movil #2
    
    df_tot <- cbind(Lag_data, ITT_data, PI_data)
    
    colnames(df_tot) <- c("Lag","ITT","PI")
    
    
    
    ###////////////////////////////////////////////////////////////////////////
    # Estimating the TOTs
    
    sum_pi_tot <- 0
    
    for (lag in df_tot$Lag+1) {
      print("Comienza lag:")
      print(lag)
      
      if (lag == 0) {
        print("Pasamos lag 0")
        
      } 
      else {
          if (lag > 2) {      
          
          pi_tot <- list() 
          pi_tot[1] <- 0
          
          for (h in 2:lag) {  
            pi_tot[h] <- df_tot$PI[h]*df_tot$TOT[lag+1-h] 
          }
          
          sum_pi_tot <- Reduce("+", pi_tot)
          print("Finaliza lag:")
          print(toString(lag))
        } 
        
        
        df_tot$TOT[lag] <-  as.numeric(df_tot$ITT[lag]) - as.numeric(sum_pi_tot)
        
         if (lag == 2) {
          pi_tot <- list()
          pi_tot[lag] <- df_tot$PI[lag]*df_tot$TOT[lag-1]
          df_tot$TOT[lag] <- as.numeric(df_tot$ITT[lag]) - as.numeric(pi_tot[lag])
          sum_pi_tot <- df_tot$TOT[lag]
          print("Finaliza lag 2")
        }
        
         if (lag == 1) {
          pi_tot <- list()
          pi_tot[lag] <- df_tot$PI[lag]*df_tot$ITT[lag]
          df_tot$TOT[lag] <- as.numeric(pi_tot[lag])
          sum_pi_tot <- df_tot$TOT[lag]
          print("Finaliza lag 1")
        }
        
      } 
    }
    
    df_tot <- within(df_tot, TOT_cumulative <- cumsum(TOT)) 
    
    ## Collecting the results data to a dataframe 
    BD_TOT_is[,ind_itt] <- df_tot$TOT   # loading individual TOTs
    BD_TOT_cs[,ind_itt] <- df_tot$TOT_cumulative  # Loading cumulative TOTs
    
    
    ind_itt <- ind_itt + 1
  } 
  
## Calculating the SE for the TOTs


## Preparing matrices to load the results. 
  BD_tot_se$Lag <- df_tot$Lag
  
  index_rest <- 1
 
 # Loading relevant data
  bd_ITTs_PIs <- read_dta(paste0(files_path,
                            "/Bootstrap/PI_",i,"_matrix_bootstrapping_LLR_500.dta"))
  bd_ITTs_PIs <- bd_ITTs_PIs[,1:num_sam+1]
  
  
## Loop to calculate the SE and collecting results
  
 while (index_rest < num_lags+1){
   
# loading the SE for cumulative TOTs
    
   # getting the average for the lag 1 of the 500 iterations
   A <- as.numeric(BD_TOT_cs[index_rest,1:num_sam])  
        # this SE comes from (https://www.stata.com/manuals/rbootstrap.pdf)
   
   # Filling the matrix of SE of the TOTc
   BD_tot_se[index_rest,3] <- sqrt(sum((A-mean(A))^2)/(length(A)-1))

# Filling the matrix with ITT SE   
 
   # getting the average for the lag 1 of the 500 iterations
   B <- as.numeric(bd_ITTs_PIs[index_rest,1:num_sam])  
        # this SE comes from (https://www.stata.com/manuals/rbootstrap.pdf)
   
   # Filling the matrix of SE of the ITTs
   BD_tot_se[index_rest,7] <- sqrt(sum((B-mean(B))^2)/(length(B)-1)) 
   
# Counting the loop
    index_rest <- index_rest+1
  }
  
  
  
###################################################################-  
#### Estimations of the cumulative TOTs ####
  
  # The code runs the estimate for each replication and then leaves a matrix of the TOTc to obtain the standard error
  
  # Loading the excel databases with the original ITT and PI estimates to calculate the original TOTc
    v_ori_PI <- read_excel(paste0(files_path,
                                   "/Full_sample/cortados_3_8/ITT_LLR_PI_cortado_3-8.xlsx"))
    v_ori_PI <- v_ori_PI[,c(2,7)]
  
  # Loading the excel databases with the original ITT and PI estimates to calculate the original TOTc
    v_ori_ITT <- read_excel(paste0(files_path,
                                   "/Full_sample/cortados_3_8/ITT_LLR_",i,"_cortado_3-8.xlsx"))
    v_ori_ITT <- v_ori_ITT[,c(2,7)]
    
  # Merging data from PI y ITT
    v_ori_PI_ITT <- merge(v_ori_PI, v_ori_ITT, 
          by = "Lag", all = TRUE)
    v_ori_PI_ITT[1,2] <- 1
    v_ori_PI_ITT[,4] <- NA
    
    names(v_ori_PI_ITT) = c("Lag","PI","ITT","TOT")
    
    head(v_ori_PI_ITT, 5)
    
    ###////////////////////////////////////////////////////////////////////////
    # Making the estimations of TOT
    
    sum_pi_tot <- 0
    
    for (lag in v_ori_PI_ITT$Lag+1) {
      print("Comienza lag:")
      print(lag)
      
      if (lag == 0) {
        print("Pasamos lag 0")
        
      } 
      else {
        if (lag > 2) {        
          
          pi_tot <- list() 
          pi_tot[1] <- 0
          
          for (h in 2:lag) {  
            pi_tot[h] <- v_ori_PI_ITT$PI[h]*v_ori_PI_ITT$TOT[lag+1-h] 
          }
          
          sum_pi_tot <- Reduce("+", pi_tot)
          print("Finaliza lag:")
          print(toString(lag))
        } 
        
        
        v_ori_PI_ITT$TOT[lag] <-  as.numeric(v_ori_PI_ITT$ITT[lag]) - as.numeric(sum_pi_tot)
      
          if (lag == 2) {
          pi_tot <- list()
          pi_tot[lag] <- v_ori_PI_ITT$PI[lag]*v_ori_PI_ITT$TOT[lag-1]
          v_ori_PI_ITT$TOT[lag] <- as.numeric(v_ori_PI_ITT$ITT[lag]) - as.numeric(pi_tot[lag])
          sum_pi_tot <- v_ori_PI_ITT$TOT[lag]
          print("Finaliza lag 2")
        }
        
         if (lag == 1) {
          pi_tot <- list()
          pi_tot[lag] <- v_ori_PI_ITT$PI[lag]*v_ori_PI_ITT$ITT[lag]
          v_ori_PI_ITT$TOT[lag] <- as.numeric(pi_tot[lag])
          sum_pi_tot <- v_ori_PI_ITT$TOT[lag]
          print("Finaliza lag 1")
        }
        
      } 
    }
  
  # Adding the cumultive TOTs variable to the base
  v_ori_PI_ITT <- within(v_ori_PI_ITT, TOT_cumulative <- cumsum(TOT))

 
  
  
  
  
    
  ###################################################################-  
  #### UNIFYING RESULTS - 95% CI ####

  # Setting the idex to 1
  index_rest <- 1

  # Adding the original TOTc and ITT data to the final results matrix and calculating the CIs
  BD_tot_se[,2] <- v_ori_PI_ITT[,5]   # TOTc
  BD_tot_se[,6] <- v_ori_PI_ITT[,3]   # ITTs
  
  # adding the 95% confidence intervals for the TOTc
  BD_tot_se[,4] <- BD_tot_se[,2] - ((1.96)*BD_tot_se[,3])
  BD_tot_se[,5] <- BD_tot_se[,2] + ((1.96)*BD_tot_se[,3])
  
  # Adding 95% confidence intervals for TOTc
  BD_tot_se[,8] <- BD_tot_se[,6] - ((1.96)*BD_tot_se[,7])
  BD_tot_se[,9] <- BD_tot_se[,6] + ((1.96)*BD_tot_se[,7])
  

  

  
  
  
  
  ###################################################################-  
  #### Graphs ITT Y TOTc ####
  
  ##/////////////////////////////////////////////////////////////////////////////////////////
  ## graph of cumulative TOTs
  plot12 <- ggplot(data=BD_tot_se,
                   aes(x=Lag, y=TOTc))+
    labs(y = paste0("TOTc ",name_y_label_list[ind_ylabel_graphs]) , x = "Lag (Months)") +
    scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96))+
    coord_cartesian(xlim= c(1.99,95))+
    theme(  text=element_text(size=12,  family="serif"),
            panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid.minor.y = element_line(size = 0.5, linetype = 'solid',
                                              colour = "azure2"),
            axis.line.x = element_line(size = 0.5, linetype = 'solid',
                                       colour = "black"),
            axis.line.y = element_line(size = 0.5, linetype = 'solid',
                                       colour = "black"),
            axis.text=element_text(size=13), 
            axis.title=element_text(size=12))+
    geom_hline(yintercept=0, linetype="dashed", color = "blue", size=0.5)+
    geom_ribbon(aes(ymin = CI_low, ymax = CI_upp), alpha = 0.2, fill = "#8B0A50")+
    geom_line(aes(Lag, CI_low), color = "grey70", size = 0.2, linetype = 2) + 
    geom_line(aes(Lag, CI_upp), color = "grey70", size = 0.2, linetype = 2) +  
    geom_line(color="grey70") +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_upp), width = 0.2, color="gray1")+ 
    geom_pointrange(data = BD_tot_se, mapping=aes(x=Lag, ymin=CI_low, ymax=CI_upp), size=0.5, color="#8B0A50") 
  plot12 
  
  ##/////////////////////////////////////////////////////////////////////////////////////////
  ## graph of ITTs  
  plot13 <- ggplot(data=BD_tot_se,
                   aes(x=Lag, y=ITT))+
    labs(y = paste0("ITT ",name_y_label_list[ind_ylabel_graphs]) , x = "Lag (Months)") +
    scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96))+
    coord_cartesian(xlim= c(1.99,95))+
    theme(  text=element_text(size=12,  family="serif"),
            panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid.minor.y = element_line(size = 0.5, linetype = 'solid',
                                              colour = "azure2"),
            axis.line.x = element_line(size = 0.5, linetype = 'solid',
                                       colour = "black"),
            axis.line.y = element_line(size = 0.5, linetype = 'solid',
                                       colour = "black"),
            axis.text=element_text(size=13), 
            axis.title=element_text(size=12))+
    geom_hline(yintercept=0, linetype="dashed", color = "blue", size=0.5)+
    geom_ribbon(aes(ymin = CI_low_itt, ymax = CI_upp_itt), alpha = 0.2, fill = "darkslategray4")+
    geom_line(aes(Lag, CI_low_itt), color = "grey70", size = 0.2, linetype = 2) + 
    geom_line(aes(Lag, CI_upp_itt), color = "grey70", size = 0.2, linetype = 2) +  
    geom_line(color="grey70") +
    geom_errorbar(aes(ymin = CI_low_itt, ymax = CI_upp_itt), width = 0.2, color="gray1")+ 
    geom_pointrange(data = BD_tot_se, mapping=aes(x=Lag, ymin=CI_low_itt, ymax=CI_upp_itt), size=0.5, color="gray1") 
  plot13
  
  # Exporting graphs
  png(file = paste0(results_path,"/Bootstrap/",i,"_Bootstrap_LLR_CI95.png"), width = 1018, height = 590)
  grid.arrange(plot13, plot12, nrow=2)
  dev.off()
  
  ind_ylabel_graphs <- ind_ylabel_graphs+1
  print(paste0("loop terminado ", i))
  
  
  ###################################################################-  
  #### Exporting tables of estimations of TOT for the appendix ####
  
  write.xlsx(BD_tot_se[,1:5], file = paste0(results_path,"/Bootstrap/",i,"TOT_AppendixTable.xlsx"), colNames = TRUE)
  
}
  
 

