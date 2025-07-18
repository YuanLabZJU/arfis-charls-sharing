setwd('~')
rm(list = ls())

library(readxl)
library(tidyverse)
library(patchwork)

# Plot the survival probability, life expectancy and life expectancy gained by Number of ARFIs

IR_adjust <- left_join(IR, Prev, by = c("age_id", "age_name")) %>% 
  left_join(ARFI_HR, by = "Number of ARFIs")

Prev_HR_Sum <- IR_adjust %>% 
  group_by(age_id) %>% 
  summarise(Prev_HR_Sum = sum(Prev_est*HR_est))

IR_adjust <- IR_adjust %>% 
  left_join(Prev_HR_Sum, by = "age_id") %>% 
  mutate(
    IR_adjust = val*HR_est/Prev_HR_Sum,
         )


life_table <- data.frame(age = rep(c(50:99),4),
                         age_id = rep(c(rep(15,5), rep(16,5), rep(17,5), rep(18,5), 
                                        rep(19,5), rep(20,5), rep(30,5), rep(31,5), 
                                        rep(32,5), rep(235,5)),4),
                         `Number of ARFIs` = c(rep(0,50),
                                               rep(1,50),
                                               rep(2,50),
                                               rep(3,50))) %>% 
  left_join(IR_adjust, by = c("age_id", "Number.of.ARFIs"="Number of ARFIs")) %>% 
  group_by(Number.of.ARFIs) %>%
  mutate(Prob_Surv = 1-IR_adjust/100000,
         Cum_Prob_Surv = cumprod(Prob_Surv),
         Life_Expectancy = sum(Cum_Prob_Surv) - cumsum(Cum_Prob_Surv)
  )

life_table$Number.of.ARFIs <- as.factor(life_table$Number.of.ARFIs)

life_table <- life_table %>% 
  group_by(age) %>% 
  mutate(Life_Expectancy_Gained = Life_Expectancy - first(Life_Expectancy))

# Survival Probability
plot_surv_prob <- ggplot(life_table, aes(x = age, y = Cum_Prob_Surv, 
                                         color = `Number.of.ARFIs`,
                                         fill = `Number.of.ARFIs`)) +
  geom_point() +
  geom_line() + 
  ggsci::scale_color_lancet() + 
  theme_bw()+labs(x = "Age (Years)", y = "Survival Probability",
                  colour = "Number of ARFIs", fill = "Number of ARFIs")

plot_le50 <- ggplot(life_table, aes(x = age, y = Life_Expectancy, 
                                    color = (`Number.of.ARFIs`),
                                    fill = (`Number.of.ARFIs`))) +
  geom_point() +
  geom_line() + 
  ggsci::scale_color_lancet() + 
  theme_bw()+labs(x = "Age (Years)", y = "Life Expectancy (Years)",
                  colour = "Number of ARFIs", fill = "Number of ARFIs")

plot_le_gained <- ggplot(life_table, aes(x = age, y = Life_Expectancy_Gained, 
                                         color = `Number.of.ARFIs`,
                                         fill = `Number.of.ARFIs`)) +
  geom_point() +
  geom_line() + 
  ggsci::scale_color_lancet() + 
  theme_bw()+labs(x = "Age (Years)", y = "Life Expectancy Loss (Years)",
                  colour = "Number of ARFIs", fill = "Number of ARFIs")

plots <- plot_surv_prob + plot_le50 + plot_le_gained + 
  plot_layout(guides = "collect",) + 
  plot_annotation( tag_level = "A") &
  theme(legend.position='bottom')

ggsave("output/LE_calc_threeversion.pdf", plots, width = 12, height = 5, units = "in", dpi = 300)




# Bootstrap to calculate CI
set.seed(123)

monte_carlo_le <- function(){
  IR_adjust <- left_join(IR, Prev, by = c("age_id", "age_name")) %>% 
    left_join(ARFI_HR, by = "Number of ARFIs") %>% 
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
      # IR_adjust = IR_mean,
    )
  life_table <- data.frame(age = rep(c(50:99),4),
                           age_id = rep(c(rep(15,5), rep(16,5), rep(17,5), rep(18,5), 
                                          rep(19,5), rep(20,5), rep(30,5), rep(31,5), 
                                          rep(32,5), rep(235,5)),4),
                           `Number of ARFIs` = c(rep(0,50),
                                                 rep(1,50),
                                                 rep(2,50),
                                                 rep(3,50))) %>% 
    left_join(IR_adjust, by = c("age_id", "Number.of.ARFIs"="Number of ARFIs")) %>% 
    group_by(Number.of.ARFIs) %>%
    mutate(Prob_Surv = 1-IR_adjust/100000,
           Cum_Prob_Surv = cumprod(Prob_Surv),
           # Life expectancy = sum of all later probabilities 
           # = sum of all Cum_Prob_Surv - sum of all previous Cum_Prob_Surv
           Life_Expectancy = sum(Cum_Prob_Surv) - cumsum(Cum_Prob_Surv)
    )
  life_table_50 <- life_table %>% 
    filter(age == 50) %>% 
    select(Number.of.ARFIs, Prob_Surv, Cum_Prob_Surv, Life_Expectancy)
  return(life_table_50)
}


# Calculate the 95% CI
apply_fun <- replicate(1000, monte_carlo_le(), simplify = FALSE)

le_50 <- bind_rows(apply_fun)
le_50_a0 <- le_50 %>% group_by(Number.of.ARFIs) %>% 
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
  select(Number.of.ARFIs, LE, YLL)
le_50_a0[1,"YLL"] <- c("0.00 (Reference)")
le_50_a0 %>% write_excel_csv("output/LE_calc_threeversion.csv")


