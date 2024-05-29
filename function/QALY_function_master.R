# Function body ----
##  scenario = "F0": 
##              The model starts with a population with NASH. 
##  scenario = "No_incidence": 
##              The model starts with a population without NASH. No risk of developing NASH (A world without NASH).
##  scenario = "No_NASH": 
##              The model starts with a population without NASH. There is a risk of developing NASH.
QALY_function_master <- function(calibrated_inc, scenario = "F0", gender = "male") {
  for (starting_age in c(5,15,25,35,45,55,65,75,90)) {
    
    age <- starting_age
    s_names <- c("No_NASH", "F0", "F1", "F2", "F3", "F4", "HCC", "DCC", "LT", "PLT", "Death")
    n_states <- length(s_names)
    n_cohort <- 1000  
    n_cycles <- 101-starting_age
    cycles <- starting_age:101
    
    # NASH incidence rate ----
    inc <- 0
    if (scenario == "No_NASH") {inc <- calibrated_inc}
    
    # Setting up the transition probability array ----
    tp.matrix<-array(data=0,dim=c(n_states,n_states,n_cycles),dimnames = list(s_names,
                                                                              s_names))
    
    for (i in 1:n_cycles) {
      
      tp.matrix["No_NASH","Death",i] <- mat_tpDn[age,gender]
      tp.matrix["No_NASH","F0",i] <- min(1,inc*incRR[age])*(1-mat_tpDn[age,gender])
      tp.matrix["No_NASH","No_NASH",i] <- (1- min(1,inc*incRR[age]))*(1 - mat_tpDn[age,gender])
      
      tp.matrix["F0","Death",i] <- mat_tpDn[age,gender]
      tp.matrix["F0","HCC",i] <- tpF0_HCC*(1-mat_tpDn[age,gender])
      tp.matrix["F0","F1",i] <- tpF0_F1*(1-mat_tpDn[age,gender])
      tp.matrix["F0","F0",i] <- (1-tpF0_HCC-tpF0_F1)*(1-mat_tpDn[age,gender])
      
      tp.matrix["F1","Death",i] <- mat_tpDn[age,gender]
      tp.matrix["F1","HCC",i] <- tpF1_HCC*(1-mat_tpDn[age,gender])
      tp.matrix["F1","F2",i] <- tpF1_F2*(1-mat_tpDn[age,gender])
      tp.matrix["F1","F1",i] <- (1-tpF1_F2-tpF1_HCC)*(1-mat_tpDn[age,gender])
      
      tp.matrix["F2","Death",i] <- mat_tpDn[age,gender]
      tp.matrix["F2","HCC",i] <- tpF2_HCC*(1-mat_tpDn[age,gender])
      tp.matrix["F2","F3",i] <- tpF2_F3*(1-mat_tpDn[age,gender])
      tp.matrix["F2","F2",i] <- (1-tpF2_HCC-tpF2_F3)*(1-mat_tpDn[age,gender])
      
      
      tp.matrix["F3","Death",i] <- mat_tpDn[age,gender]
      tp.matrix["F3","HCC",i] <- tpF3_HCC*(1-mat_tpDn[age,gender])
      tp.matrix["F3","F4",i] <- tpF3_F4*(1-mat_tpDn[age,gender])
      tp.matrix["F3","F3",i] <- (1-tpF3_HCC-tpF3_F4)*(1-mat_tpDn[age,gender])
      
      tp.matrix["F4","Death",i] <- min(mat_tpDn[age,gender]+tpF4_Death,1) 
      tp.matrix["F4","DCC",i] <- tpF4_DCC*(1- min(mat_tpDn[age,gender]+tpF4_Death,1))
      tp.matrix["F4","HCC",i] <- tpF4_HCC*(1- min(mat_tpDn[age,gender]+tpF4_Death,1))
      tp.matrix["F4","F4",i] <- (1-tpF4_DCC-tpF4_HCC)*(1- min(mat_tpDn[age,gender]+tpF4_Death,1))
      
      tp.matrix["DCC","Death",i] <- min(mat_tpDn[age,gender]+tpDCC_Death,1)
      tp.matrix["DCC","LT",i] <- tpDCC_LT*(1-min(mat_tpDn[age,gender]+tpDCC_Death,1))
      tp.matrix["DCC","HCC",i] <- tpDCC_HCC*(1-min(mat_tpDn[age,gender]+tpDCC_Death,1))
      tp.matrix["DCC","DCC",i] <- (1-tpDCC_LT-tpDCC_HCC)*(1-min(mat_tpDn[age,gender]+tpDCC_Death,1))
      
      tp.matrix["HCC","Death",i] <- min(mat_tpDn[age,gender]+tpHCC_Death,1)
      tp.matrix["HCC","LT",i] <-  (tpHCC_LT)*(1-min(mat_tpDn[age,gender]+tpHCC_Death,1))
      tp.matrix["HCC","HCC",i] <- (1-tpHCC_LT)*(1- min(mat_tpDn[age,gender]+tpHCC_Death,1))
      
      tp.matrix["LT","Death",i] <- min(mat_tpDn[age,gender]+tpLT_Death,1)
      tp.matrix["LT","PLT",i] <- (1 - min(mat_tpDn[age,gender]+tpLT_Death,1))
      
      tp.matrix["PLT","Death",i] <- min(mat_tpDn[age,gender]+tpPLT_Death,1)
      tp.matrix["PLT","PLT",i] <- (1 - min(mat_tpDn[age,gender]+tpPLT_Death,1))
      
      tp.matrix["Death","Death",i] <- 1
      
      
      age <- min(age+1,101)
    }
    
    
    pop <- array(data = NA,
                 dim = c(n_cycles, n_states),
                 dimnames = list(NULL, state = s_names))
    
    # Seeding the starting state ----
    pop[1,"No_NASH"] <- switch(scenario,
                               "F0" = 0,
                               "No_incidence" = n_cohort,
                               "No_NASH" = n_cohort)
    pop[1,"F0"] <- switch(scenario,
                          "F0" = n_cohort,
                          "No_incidence" = 0,
                          "No_NASH" = 0)
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
    for (i in 2:n_cycles) {
      pop[i,] <- pop[i-1,]%*%tp.matrix[,,i-1]
    }
    
    pop_t <- t(pop) # transposing the matrix trace so that it can be multiplied by c(1,1,1,1,1,1,1,1,1,0) to calculate LE
    # pop_t
    # dim(pop_t)
    
    # cycle_empty_array <-
    #   array(NA,
    #         dim = c(1, n_cycles),
    #         dimnames = list(1,
    #                         cycle = NULL))       # setting up matrix
    # 
    # 
    # LE <- cycle_empty_array
    
    LE <- paste("LE", starting_age, sep = "")
    # multiply matrix of 1s for non-death states & 0 for death state by matrix trace
    assign(LE, c(1,1,1,1,1,1,1,1,1,1,0) %*% pop_t[,] )
    
    
    # QALY ----
    
    # setting up matrix
    q_matrix <- array(data = NA,
                      dim = c(n_cycles, n_states),
                      dimnames = list(NULL, state = s_names))
    
    q_matrix_2 <- array(data = NA,
                        dim = c(n_cycles, n_states),
                        dimnames = list(NULL, state = s_names))
    
    
    # creating the matrix of utilities which vary by age 
    age <- starting_age
    for (i in 1:n_cycles) {
      
      q_matrix[i,"No_NASH"] <- mat_utility[age,gender]
      q_matrix[i,"F0"] <- mat_utility[age,gender]
      q_matrix[i,"F1"] <- mat_utility[age,gender]
      q_matrix[i,"F2"] <- mat_utility[age,gender]
      q_matrix[i,"F3"] <- mat_utility[age,gender]
      q_matrix[i,"F4"] <- mat_utility[age,gender]
      q_matrix[i,"HCC"] <- mat_utility[age,gender]
      q_matrix[i,"DCC"] <-  mat_utility[age,gender]
      q_matrix[i,"LT"] <- mat_utility[age,gender]
      q_matrix[i,"PLT"] <- mat_utility[age,gender]
      q_matrix[i,"Death"] <- 0
      
      age <- min(age+1,101)
    }
    
    
    # utility ratios
    age <-starting_age
    for (i in 1:n_cycles) {
      
      q_matrix_2[i,"No_NASH"] <- 1
      q_matrix_2[i,"F0"] <- UD_F0
      q_matrix_2[i,"F1"] <- UD_F1
      q_matrix_2[i,"F2"] <- UD_F2
      q_matrix_2[i,"F3"] <- UD_F3
      q_matrix_2[i,"F4"] <- UD_F4
      q_matrix_2[i,"HCC"] <- UD_HCC
      q_matrix_2[i,"DCC"] <-  UD_DCC
      q_matrix_2[i,"LT"] <- UD_LT
      q_matrix_2[i,"PLT"] <- UD_PLT
      q_matrix_2[i,"Death"] <- 0
      
      age <- min(age+1,101)
    }
    
    # cross multiply matrices
    q_matrix_3 <-    q_matrix*q_matrix_2                    
    
    
    # quality-adjusted life expectancy with utility decrements and no discounting
    
    cycle_QALE <- paste("cycle_QALE", starting_age, sep = "")
    assign(cycle_QALE, q_matrix_3*pop[, ]  )
    # cycle_QALE<-    q_matrix_3*pop[, ]                       #multiplying the qalys by the matrix trace 
    # cycle_QALE
    
  }
  
  
  # Calculating total LE
  
  total_LE5 <- sum(LE5[1, -1])            
  total_LE15 <- sum(LE15[1, -1])            
  total_LE25 <- sum(LE25[1, -1])            
  total_LE35 <- sum(LE35[1, -1])            
  total_LE45 <- sum(LE45[1, -1])            
  total_LE55 <- sum(LE55[1, -1])            
  total_LE65 <- sum(LE65[1, -1])            
  total_LE75 <- sum(LE75[1, -1])            
  total_LE90 <- sum(LE90[1, -1])            
  
  
  # Calculating total QALYs
  # summing columns to get total QALYs by state for each cycle 
  
  cycle_QALYs5 <- colSums(cycle_QALE5)        
  cycle_QALYs15 <- colSums(cycle_QALE15)
  cycle_QALYs25 <- colSums(cycle_QALE25)
  cycle_QALYs35 <- colSums(cycle_QALE35)
  cycle_QALYs45 <- colSums(cycle_QALE45)
  cycle_QALYs55 <- colSums(cycle_QALE55)
  cycle_QALYs65 <- colSums(cycle_QALE65)
  cycle_QALYs75 <- colSums(cycle_QALE75)
  cycle_QALYs90 <- colSums(cycle_QALE90)
  
  # summing total QALYs
  total_QALYs5<- sum(cycle_QALYs5)             
  total_QALYs15<- sum(cycle_QALYs15)             
  total_QALYs25<- sum(cycle_QALYs25)             
  total_QALYs35<- sum(cycle_QALYs35)             
  total_QALYs45<- sum(cycle_QALYs45)             
  total_QALYs55<- sum(cycle_QALYs55)             
  total_QALYs65<- sum(cycle_QALYs65)             
  total_QALYs75<- sum(cycle_QALYs75)             
  total_QALYs90<- sum(cycle_QALYs90)  
  
  # print(total_QALYs5)
  # print(total_QALYs15)
  # print(total_QALYs25)
  # print(total_QALYs35)
  # print(total_QALYs45)
  # print(total_QALYs55)
  # print(total_QALYs65)
  # print(total_QALYs75)
  # print(total_QALYs90)
  # 
  # print(total_LE5)
  # print(total_LE15)
  # print(total_LE25)
  # print(total_LE35)
  # print(total_LE45)
  # print(total_LE55)
  # print(total_LE65)
  # print(total_LE75)
  # print(total_LE90)
  
  
  QALYbyAge <- c(total_QALYs5, 
                 total_QALYs15, 
                 total_QALYs25,             
                 total_QALYs35,             
                 total_QALYs45,             
                 total_QALYs55,             
                 total_QALYs65,             
                 total_QALYs75,             
                 total_QALYs90)
  
  LEbyAge <- c(total_LE5, 
               total_LE15, 
               total_LE25,             
               total_LE35,             
               total_LE45,             
               total_LE55,             
               total_LE65,             
               total_LE75,             
               total_LE90)
  
  age_range <- c("0-9","10-19","20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  return.mat <- data.frame(age_range,QALYbyAge,LEbyAge)
  
  
  ##################################################
  
  # Calculating lifetime burden by starting age
  i <- 8
  
  for (starting_age in c(5,15,25,35,45,55,65,75,90)) {
    age <- 1:(101-starting_age)
    age_grp <- cut(age, breaks=c(0,seq(from = 9, by = 10, length.out = i), (101-starting_age)))
    
    # QALY by age group
    
    cycle_QALE_agegrp <- paste("cycle_QALE_agegrp", starting_age, sep = "")
    cycle_QALE <- paste("cycle_QALE", starting_age, sep = "")
    assign(cycle_QALE_agegrp, rowsum(get(cycle_QALE), age_grp))
    
    cycle_QALY_agegrp <- paste("cycle_QALY_agegrp", starting_age, sep = "")
    assign(cycle_QALY_agegrp, rowSums(get(cycle_QALE_agegrp)))
    
    temp <- data.matrix(get(cycle_QALY_agegrp))
    temp <- rbind(matrix(rep(0, times = 8-i), ncol = 1),temp)
    rownames(temp) <- age_range 
    
    cycle_QALY <- paste("cycle_QALY", starting_age, sep = "")
    assign(cycle_QALY, temp)
    
    # print(cycle_QALY)
    # print(get(cycle_QALY))
    
    # LE by age group
    age_grp <- cut(age, breaks=c(1,seq(from = 9, by = 10, length.out = i), (101-starting_age)))
    age_grp <- na.omit(age_grp)
    
    LE_t <- paste("LE_t", starting_age, sep = "")
    LE <- paste("LE", starting_age, sep = "")
    temp <- t(get(LE))
    temp <- temp[-1,]
    assign(LE_t, temp)
    
    temp <- rowsum(temp, age_grp)
    temp <- rbind(matrix(rep(0, times = 8-i), ncol = 1),temp)
    rownames(temp) <- age_range
    cycle_LE <- paste("cycle_LE", starting_age, sep = "")
    assign(cycle_LE, temp)
    
    # print(cycle_LE)
    # print(get(cycle_LE))
    
    return.cycle <- paste("return.cycle", starting_age, sep = "")
    assign(return.cycle, data.frame(rownames(get(cycle_QALY)),get(cycle_QALY),get(cycle_LE)))
    
    i <- i - 1
  }
  
  # Population US ----
  
  cycle_LE_0_9 <- cycle_LE5["0-9",1]
  cycle_LE_10_19 <- cycle_LE5["10-19",1] + cycle_LE15["10-19",1]
  cycle_LE_20_29 <- cycle_LE5["20-29",1] + cycle_LE15["20-29",1] + cycle_LE25["20-29",1]
  cycle_LE_30_39 <- cycle_LE5["30-39",1] + cycle_LE15["30-39",1] + cycle_LE25["30-39",1] + cycle_LE35["30-39",1]
  cycle_LE_40_49 <- cycle_LE5["40-49",1] + cycle_LE15["40-49",1] + cycle_LE25["40-49",1] + cycle_LE35["40-49",1] + cycle_LE45["40-49",1]
  cycle_LE_50_59 <- cycle_LE5["50-59",1] + cycle_LE15["50-59",1] + cycle_LE25["50-59",1] + cycle_LE35["50-59",1] + cycle_LE45["50-59",1] + cycle_LE55["50-59",1]
  cycle_LE_60_69 <- cycle_LE5["60-69",1] + cycle_LE15["60-69",1] + cycle_LE25["60-69",1] + cycle_LE35["60-69",1] + cycle_LE45["60-69",1] + cycle_LE55["60-69",1] + cycle_LE65["60-69",1]
  cycle_LE_70_79 <- cycle_LE5["70-79",1] + cycle_LE15["70-79",1] + cycle_LE25["70-79",1] + cycle_LE35["70-79",1] + cycle_LE45["70-79",1] + cycle_LE55["70-79",1] + cycle_LE65["70-79",1] + cycle_LE75["70-79",1]
  cycle_LE_80plus <- cycle_LE5["80+",1] + cycle_LE15["80+",1] + cycle_LE25["80+",1] + cycle_LE35["80+",1] + cycle_LE45["80+",1] + cycle_LE55["80+",1] + cycle_LE65["80+",1] + cycle_LE75["80+",1] + cycle_LE90["80+",1]
  
  
  cycle_QALY_0_9 <- cycle_QALY5["0-9",1]
  cycle_QALY_10_19 <- cycle_QALY5["10-19",1] + cycle_QALY15["10-19",1]
  cycle_QALY_20_29 <- cycle_QALY5["20-29",1] + cycle_QALY15["20-29",1] + cycle_QALY25["20-29",1]
  cycle_QALY_30_39 <- cycle_QALY5["30-39",1] + cycle_QALY15["30-39",1] + cycle_QALY25["30-39",1] + cycle_QALY35["30-39",1]
  cycle_QALY_40_49 <- cycle_QALY5["40-49",1] + cycle_QALY15["40-49",1] + cycle_QALY25["40-49",1] + cycle_QALY35["40-49",1] + cycle_QALY45["40-49",1]
  cycle_QALY_50_59 <- cycle_QALY5["50-59",1] + cycle_QALY15["50-59",1] + cycle_QALY25["50-59",1] + cycle_QALY35["50-59",1] + cycle_QALY45["50-59",1] + cycle_QALY55["50-59",1]
  cycle_QALY_60_69 <- cycle_QALY5["60-69",1] + cycle_QALY15["60-69",1] + cycle_QALY25["60-69",1] + cycle_QALY35["60-69",1] + cycle_QALY45["60-69",1] + cycle_QALY55["60-69",1] + cycle_QALY65["60-69",1]
  cycle_QALY_70_79 <- cycle_QALY5["70-79",1] + cycle_QALY15["70-79",1] + cycle_QALY25["70-79",1] + cycle_QALY35["70-79",1] + cycle_QALY45["70-79",1] + cycle_QALY55["70-79",1] + cycle_QALY65["70-79",1] + cycle_QALY75["70-79",1]
  cycle_QALY_80plus <- cycle_QALY5["80+",1] + cycle_QALY15["80+",1] + cycle_QALY25["80+",1] + cycle_QALY35["80+",1] + cycle_QALY45["80+",1] + cycle_QALY55["80+",1] + cycle_QALY65["80+",1] + cycle_QALY75["80+",1] + cycle_QALY90["80+",1]
  
  # US male population
  pop_gender <- setNames(mat_uspop_10[,gender],mat_uspop_10[,"agegroup"])
  
  # calculated LEs and QALYs multiplied by the US population in each age group (I divided by 1000 because the initial cohort size when running the model was set to 1000)
  
  pop_LE_0_9  <- (cycle_LE_0_9/1000)*pop_gender["0-9"] 
  pop_LE_10_19  <- (cycle_LE_10_19/1000)*pop_gender["10-19"] 
  pop_LE_20_29  <- (cycle_LE_20_29/1000)*pop_gender["20-29"] 
  pop_LE_30_39  <- (cycle_LE_30_39/1000)*pop_gender["30-39"] 
  pop_LE_40_49  <- (cycle_LE_40_49/1000)*pop_gender["40-49"] 
  pop_LE_50_59  <- (cycle_LE_50_59/1000)*pop_gender["50-59"] 
  pop_LE_60_69  <- (cycle_LE_60_69/1000)*pop_gender["60-69"] 
  pop_LE_70_79  <- (cycle_LE_70_79/1000)*pop_gender["70-79"] 
  pop_LE_80plus  <- (cycle_LE_80plus/1000)*pop_gender["80+"] 
  
  
  pop_QALY_0_9  <- (cycle_QALY_0_9/1000)*pop_gender["0-9"] 
  pop_QALY_10_19  <- (cycle_QALY_10_19/1000)*pop_gender["10-19"] 
  pop_QALY_20_29  <- (cycle_QALY_20_29/1000)*pop_gender["20-29"] 
  pop_QALY_30_39  <- (cycle_QALY_30_39/1000)*pop_gender["30-39"] 
  pop_QALY_40_49  <- (cycle_QALY_40_49/1000)*pop_gender["40-49"] 
  pop_QALY_50_59  <- (cycle_QALY_50_59/1000)*pop_gender["50-59"] 
  pop_QALY_60_69  <- (cycle_QALY_60_69/1000)*pop_gender["60-69"] 
  pop_QALY_70_79  <- (cycle_QALY_70_79/1000)*pop_gender["70-79"] 
  pop_QALY_80plus  <- (cycle_QALY_80plus/1000)*pop_gender["80+"] 
  
  
  
  
  ######################
  
  
  QALY_pop <- c(pop_QALY_0_9, 
                pop_QALY_10_19, 
                pop_QALY_20_29, 
                pop_QALY_30_39, 
                pop_QALY_40_49, 
                pop_QALY_50_59, 
                pop_QALY_60_69, 
                pop_QALY_70_79, 
                pop_QALY_80plus)
  
  LE_pop <- c(pop_LE_0_9, 
              pop_LE_10_19, 
              pop_LE_20_29, 
              pop_LE_30_39, 
              pop_LE_40_49, 
              pop_LE_50_59, 
              pop_LE_60_69, 
              pop_LE_70_79, 
              pop_LE_80plus)
  
  return.pop <- data.frame(age_range,QALY_pop,LE_pop)
  
  
  # # function will print the estimates   
  # print(pop_QALY_0_9)
  # print(pop_QALY_10_19)
  # print(pop_QALY_20_29)
  # print(pop_QALY_30_39)
  # print(pop_QALY_40_49)
  # print(pop_QALY_50_59)
  # print(pop_QALY_60_69)
  # print(pop_QALY_70_79)
  # print(pop_QALY_80plus)
  
  # print(pop_LE_0_9)
  # print(pop_LE_10_19)
  # print(pop_LE_20_29)
  # print(pop_LE_30_39)
  # print(pop_LE_40_49)
  # print(pop_LE_50_59)
  # print(pop_LE_60_69)
  # print(pop_LE_70_79)
  # print(pop_LE_80plus)
  
  
  # Return a list ----
  return(list(return.mat = return.mat, 
              
              return.cycle5 = return.cycle5,
              return.cycle15 = return.cycle15,
              return.cycle25 = return.cycle25,
              return.cycle35 = return.cycle35,
              return.cycle45 = return.cycle45,
              return.cycle55 = return.cycle55,
              return.cycle65 = return.cycle65,
              return.cycle75 = return.cycle75,
              return.cycle90 = return.cycle90, 
              
              return.pop = return.pop
  )
  )
}
