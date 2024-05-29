##

CalibrateInc2 <- function(inc,
                          targetpre,
                          F0 = FALSE, F1 = FALSE, F2 = FALSE, F3 = FALSE, F4 = FALSE, HCC = FALSE, DCC = FALSE, LT = FALSE, PLT = FALSE, 
                          AdultOnly = FALSE, 
                          ExactNum = FALSE){
  

  

  
# F0 = F
# F1 = F
# F2 = F
# F3 = F
# F4 = F
# HCC = F
# DCC = F
# LT = F
# PLT = F
# AdultOnly = F
# ExactNum = FALSE
# 
  # F0 = F
  # F1 = F
  # F2 = F
  # F3 = F
  # F4 = T
  # HCC = F
  # DCC = T
  # LT = F
  # PLT = F
  # AdultOnly = T
  # ExactNum = F
  
  # # The parameter to be calibrated: NASH incidence rate ----
  # inc <- 0.000452466
  
  s_names <- c("No_NASH", "F0", "F1", "F2", "F3", "F4", "HCC", "DCC", "LT", "PLT", "Death")
  n_states <- length(s_names)
  agegroups <- seq(from = 0, to = 85, by = 5)
  n_agegroups <- length(agegroups)
  
  NoP_calculated <- array(data = 0, dim = c(n_agegroups,n_states - 2), dimnames = list(agegroup = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+"),
                                                                                       s_names[2:(n_states-1)]))
  
  StatesIncl <- as.numeric(c(F0, F1, F2, F3, F4, HCC, DCC, LT, PLT))
  

  ### For current living cohort ----
  for (gender in c("male","female")) {
    # gender <- "male"
    
    
    counter <- 0
    for (starting_age in agegroups) {
      # starting_age = 0
      counter <- counter + 1
      age <- starting_age
      
      n_cohort <- mat_uspop[counter,gender] 
      n_cycles <- 88-starting_age
      
      # number of newly developed F0 patients
      NoP_calculated[counter,"F0"] <- NoP_calculated[counter,"F0"] + n_cohort * min(1,inc*incRR[max(age,1)])
      
      # number of survived patients
      if (counter < n_agegroups) {
        # Setting up the transition probability array ----
        tp.matrix<-array(data=0,dim=c(n_states,n_states,n_cycles),dimnames = list(s_names,
                                                                                  s_names))
        
        for (i in 1:n_cycles) {
          
          tp.matrix["No_NASH","Death",i] <- mat_tpDn[max(age,1),gender]
          tp.matrix["No_NASH","F0",i] <- min(1,inc*incRR[max(age,1)])*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["No_NASH","No_NASH",i] <- (1- min(1,inc*incRR[max(age,1)]))*(1 - mat_tpDn[max(age,1),gender])
          
          tp.matrix["F0","Death",i] <- mat_tpDn[max(age,1),gender]
          tp.matrix["F0","HCC",i] <- tpF0_HCC*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F0","F1",i] <- tpF0_F1*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F0","F0",i] <- (1-tpF0_HCC-tpF0_F1)*(1-mat_tpDn[max(age,1),gender])
          
          tp.matrix["F1","Death",i] <- mat_tpDn[max(age,1),gender]
          tp.matrix["F1","HCC",i] <- tpF1_HCC*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F1","F2",i] <- tpF1_F2*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F1","F1",i] <- (1-tpF1_F2-tpF1_HCC)*(1-mat_tpDn[max(age,1),gender])
          
          tp.matrix["F2","Death",i] <- mat_tpDn[max(age,1),gender]
          tp.matrix["F2","HCC",i] <- tpF2_HCC*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F2","F3",i] <- tpF2_F3*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F2","F2",i] <- (1-tpF2_HCC-tpF2_F3)*(1-mat_tpDn[max(age,1),gender])
          
          
          tp.matrix["F3","Death",i] <- mat_tpDn[max(age,1),gender]
          tp.matrix["F3","HCC",i] <- tpF3_HCC*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F3","F4",i] <- tpF3_F4*(1-mat_tpDn[max(age,1),gender])
          tp.matrix["F3","F3",i] <- (1-tpF3_HCC-tpF3_F4)*(1-mat_tpDn[max(age,1),gender])
          
          tp.matrix["F4","Death",i] <- min(mat_tpDn[max(age,1),gender]+tpF4_Death,1) 
          tp.matrix["F4","DCC",i] <- tpF4_DCC*(1- min(mat_tpDn[max(age,1),gender]+tpF4_Death,1))
          tp.matrix["F4","HCC",i] <- tpF4_HCC*(1- min(mat_tpDn[max(age,1),gender]+tpF4_Death,1))
          tp.matrix["F4","F4",i] <- (1-tpF4_DCC-tpF4_HCC)*(1- min(mat_tpDn[max(age,1),gender]+tpF4_Death,1))
          
          tp.matrix["DCC","Death",i] <- min(mat_tpDn[max(age,1),gender]+tpDCC_Death,1)
          tp.matrix["DCC","LT",i] <- tpDCC_LT*(1-min(mat_tpDn[max(age,1),gender]+tpDCC_Death,1))
          tp.matrix["DCC","HCC",i] <- tpDCC_HCC*(1-min(mat_tpDn[max(age,1),gender]+tpDCC_Death,1))
          tp.matrix["DCC","DCC",i] <- (1-tpDCC_LT-tpDCC_HCC)*(1-min(mat_tpDn[max(age,1),gender]+tpDCC_Death,1))
          
          tp.matrix["HCC","Death",i] <- min(mat_tpDn[max(age,1),gender]+tpHCC_Death,1)
          tp.matrix["HCC","LT",i] <-  (tpHCC_LT)*(1-min(mat_tpDn[max(age,1),gender]+tpHCC_Death,1))
          tp.matrix["HCC","HCC",i] <- (1-tpHCC_LT)*(1- min(mat_tpDn[max(age,1),gender]+tpHCC_Death,1))
          
          tp.matrix["LT","Death",i] <- min(mat_tpDn[max(age,1),gender]+tpLT_Death,1)
          tp.matrix["LT","PLT",i] <- (1 - min(mat_tpDn[max(age,1),gender]+tpLT_Death,1))
          
          tp.matrix["PLT","Death",i] <- min(mat_tpDn[max(age,1),gender]+tpPLT_Death,1)
          tp.matrix["PLT","PLT",i] <- (1 - min(mat_tpDn[max(age,1),gender]+tpPLT_Death,1))
          
          tp.matrix["Death","Death",i] <- 1
          
          
          age <- min(age+1,101)
        }
        
        
        pop <- array(data = NA,
                     dim = c(n_cycles, n_states+2),
                     dimnames = list(NULL, state = c("Year","Age",s_names)))
        
        # Seeding the starting state ----
        pop[1,"Year"] <- 0
        pop[1,"Age"] <- starting_age
        
        pop[1,"No_NASH"] <- n_cohort
        pop[1,"F0"] <- 0
        pop[1,"F1"] <- 0
        pop[1,"F2"] <- 0
        pop[1,"F3"] <- 0
        pop[1,"F4"] <- 0
        pop[1,"HCC"] <- 0
        pop[1,"DCC"] <- 0
        pop[1,"LT"] <- 0
        pop[1,"PLT"] <- 0
        pop[1,"Death"] <- 0
        
        # Population matrix trace ----
        for (i in 2:5) {
          pop[i,1] <- pop[i-1,1] + 1
          pop[i,2] <- pop[i-1,2] + 1
          pop[i,-(1:2)] <- pop[i-1,-(1:2)]%*%tp.matrix[,,i-1]
        }
        

        for (i in 6:n_cycles) {
          pop[i,1] <- pop[i-1,1] + 1
          pop[i,2] <- pop[i-1,2] + 1
          pop[i,-(1:3)] <- pop[i-1,-(1:3)]%*%tp.matrix[-1,-1,i-1]
        }
        
        for (i in (counter+1):n_agegroups) {
          for (s in s_names[2:10]) {
            NoP_calculated[i,s] <- NoP_calculated[i,s] + pop[pop[,"Age"]==5*i-3,s]
            
            
          }
        }
      }
      
    }
    
  }
  
  # return the prevalence
  State_NoP_calculated <- NoP_calculated %*% StatesIncl
  
  SUMNoP_calculated <- ifelse(AdultOnly,0.4*State_NoP_calculated[4]+sum(State_NoP_calculated[-(1:4)]),sum(State_NoP_calculated))
  Prev_calculated <- SUMNoP_calculated / ifelse(ExactNum,1,ifelse(AdultOnly,0.4*sum(mat_uspop[4,2:3])+sum(mat_uspop[-(1:4),2:3]),sum(mat_uspop[,2:3])))
  
  prev_byage <- State_NoP_calculated / rowSums(mat_uspop[,2:3])
  
  TBmin <- (Prev_calculated - targetpre)^2
  return(TBmin)
  # return(Prev_calculated = Prev_calculated)
}


