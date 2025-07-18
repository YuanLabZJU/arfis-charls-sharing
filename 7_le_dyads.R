setwd('D:/Dataset/CHARLS')
rm(list = ls())

library(tidyverse)
library(tableone)
library(survival)
library(readxl)

load('data/alldata_arfi_imputation_Vad.rda')

arfis <- c('r1vi', 'r1hi', 'r1ci', 'r1sd', 'r1dp', 'r1ad')
IR <- read_excel(paste0('data/Life Expectancy - CI.xlsx'), sheet = "MortalityRates")

monte_carlo_le <- function() {
  IR_adjust <- left_join(IR, Prev, by = c("age_id", "age_name")) %>% 
    left_join(ARFI_HR, by = "Exposure") %>% 
    rowwise() %>% 
    mutate(IR_mean = val,
           IR_sd = (val-lower)/1.96,
           IR_rand = rnorm(1,IR_mean,IR_sd),
           
           Prev_mean = Prev_est,
           Prev_sd = (Prev_est-Prev_low)/1.96,
           Prev_rand = rnorm(1,Prev_mean,Prev_sd),
           
           logHR_mean = log(HR_est),
           logHR_sd = (log(HR_est)-log(HR_low))/1.96,
           logHR_rand = rnorm(1,logHR_mean,logHR_sd),
           HR_rand = exp(logHR_rand)
    )
  Prev_HR_Sum <- IR_adjust %>% 
    group_by(age_id) %>% 
    summarise(Prev_HR_Sum = sum(Prev_rand*HR_rand))
  IR_adjust <- IR_adjust %>% 
    left_join(Prev_HR_Sum, by = "age_id") %>% 
    mutate(
      IR_adjust = IR_mean*HR_rand/Prev_HR_Sum,
    )
  life_table <- data.frame(age = rep(c(50:99),3),
                           age_id = rep(c(rep(15,5), rep(16,5), rep(17,5), rep(18,5), 
                                          rep(19,5), rep(20,5), rep(30,5), rep(31,5), 
                                          rep(32,5), rep(235,5)), 3),
                           Exposure = c(rep(0,50), rep(1,50), rep(2,50))) %>% 
    left_join(IR_adjust, by = c("age_id", "Exposure")) %>% 
    group_by(Exposure) %>%
    mutate(Prob_Surv = 1-IR_adjust/100000,
           Cum_Prob_Surv = cumprod(Prob_Surv),
           Life_Expectancy = sum(Cum_Prob_Surv) - cumsum(Cum_Prob_Surv)
    )
  life_table_50 <- life_table %>% 
    filter(age == 50) %>% 
    select(Exposure, Prob_Surv, Cum_Prob_Surv, Life_Expectancy)
  return(life_table_50)
}

le_dyads <- c('', '', '')

for (i in 1:5) {
  
  for (j in (i+1):6){
    
    message(arfis[i])
    message(arfis[j])
    pairdata <- alldata %>% 
      mutate(
        e1 = alldata[[arfis[i]]], 
        e2 = alldata[[arfis[j]]], 
        exposure = e1+e2, 
        OtherARFIs = rowSums(alldata[arfis[-c(i,j)]], na.rm = T), 
        OtherARFIs = ifelse(OtherARFIs>0, 1, OtherARFIs)
      )
    
    pre_data <- pairdata %>% 
      filter(!is.na(exposure), !is.na(mo_oc)) %>% 
      group_by(exposure, AGC) %>% 
      summarise(
        n = n(), 
        pre = sum(weight*mo_oc) / sum(weight), 
        lower = pre-sqrt(pre*(1-pre)/n), 
        upper = pre+sqrt(pre*(1-pre)/n)
      )
    
    pair_cox <- coxph(Surv(mo_t2-2011, mo_oc)~factor(exposure)+
                        AGE+GENDER+EDU+MARR+HUKOU+SET+INCOME+BMI+SMK+DRK+STR+CAN+HBP+DIA+HEA+OtherARFIs, 
                      data=pairdata)
    summary_results <- ShowRegTable(pair_cox, printToggle = F)
    
    ARFI_HR <- data.frame(
      Exposure = c(0, 1, 2), 
      HR_est = c(1, substr(summary_results[1, 1], 1, 4), substr(summary_results[2, 1], 1, 4)), 
      HR_low = c(1, substr(summary_results[1, 1], 7, 10), substr(summary_results[2, 1], 7, 10)), 
      HR_hi = c(1, substr(summary_results[1, 1], 13, 16), substr(summary_results[2, 1], 13, 16))
    ) %>% 
      mutate(
        across(HR_est:HR_hi, as.numeric)
      )
    
    Prev <- data.frame(
      age_id = rep(IR$age_id, 3), 
      age_name = rep(IR$age_name, 3), 
      Exposure = rep(0:2, each = 10), 
      Prev_est = pre_data[c(1:8, 8, 8, 9:16, 16, 16, 17:24, 24, 24), ]$pre, 
      Prev_low = pre_data[c(1:8, 8, 8, 9:16, 16, 16, 17:24, 24, 24), ]$lower, 
      Prev_hi = pre_data[c(1:8, 8, 8, 9:16, 16, 16, 17:24, 24, 24), ]$upper
    )
    
    
    set.seed(42)
    
    message('Monte Carlo begins...')
    apply_fun <- replicate(1000, monte_carlo_le(), simplify = FALSE)
    message('Monte Carlo finished!')
    
    le_50 <- bind_rows(apply_fun)
    le_50_a0 <- le_50 %>% group_by(Exposure) %>% 
      summarise(LE_mean = mean(Life_Expectancy),
                LE_sd = sd(Life_Expectancy),
                LE_lower = LE_mean-qt(0.975,999)*LE_sd,
                LE_upper = LE_mean+qt(0.975,999)*LE_sd,
                LE = paste0(round(LE_mean,2), " (", round(LE_lower,2), ", ", 
                            round(LE_upper,2), ")")) %>%
      ungroup() %>% 
      mutate(
        YLL_mean = LE_mean-first(LE_mean),
        YLL_se = sqrt(LE_sd^2 + first(LE_sd)^2),
        YLL_lower = YLL_mean-qt(0.975,999)*YLL_se,
        YLL_upper = YLL_mean+qt(0.975,999)*YLL_se,
        YLL = paste0(round(YLL_mean,2), " (", round(YLL_lower,2), ", ", 
                     round(YLL_upper,2), ")")) %>% 
      select(Exposure, LE, YLL)
    
    le_50_a0[1,"YLL"] <- c("0.00 (Reference)")
    print(le_50_a0)
    
    le_dyads <- rbind(
      le_dyads, 
      c(arfis[i], arfis[j], ''),
      le_50_a0
    )
    
  }
  
}

write.csv(le_50_a0, file = 'data/le_dyads.csv')







