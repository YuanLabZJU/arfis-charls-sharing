setwd('~')
rm(list = ls())
gc()

library(tidyverse)
library(segmented)
library(readxl)
library(patchwork)
load("data/seg-prov.rda")

CHARLS_YB_Census_2020 <- read_excel("CHARLS-YB-Census-2020.xlsx")

load("2020byprovince.rda")


weight_rep <- read_excel("Weight-Rep-2011.xlsx")
names(weight_rep)
ggplot(weight_rep, aes(x = `% of National 50+ Population (6th National Census)`,
                       y = `% of weight (study population)`)) + 
  geom_text(aes(label = `Province-level region`)) + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw()

sums <- alldata %>% 
  group_by(province) %>%
  summarise(
    n = n(),
    wgts = sum(weight, na.rm = TRUE),
    
    vi_prev = weighted.mean(r5vi, weight, na.rm = TRUE),
    hi_prev = weighted.mean(r5hi, weight, na.rm = TRUE),
    sd_prev = weighted.mean(r5sd, weight, na.rm = TRUE),
    ci_prev = weighted.mean(r5ci, weight, na.rm = TRUE),
    dp_prev = weighted.mean(r5dp, weight, na.rm = TRUE),
    ad_prev = weighted.mean(r5ad, weight, na.rm = TRUE),
  ) %>% 
  ungroup() %>% 
  left_join(CHARLS_YB_Census_2020, by = c("province" = "CHARLS Name")) %>% 
  mutate(log_gdppc = log10(`人均GDP`),
         wgts_rel = wgts/sum(wgts)*sum(n),
         pop50 = `人口统计-50岁以上人口_百分比`
  ) %>% 
  filter(n>=200)

plot_break <- function(outcome_name){
  outcome<- pull(sums, outcome_name)
  glm_fi <- glm(outcome~log_gdppc + pop50 + I(pop50^2), data=sums, 
                family = binomial(link = "log")
                , weights = wgts_rel
  )
  summary(glm_fi)
  
  seg_fi <- segmented(glm_fi, seg.Z = ~ log_gdppc)
  sum_seg <- summary(seg_fi)
  davies <- davies.test(glm_fi, seg.Z = ~ log_gdppc)
  break_pt <- davies[["statistic"]][["'best' at"]]
  p_break <- davies[["p.value"]]
  p_op <- ifelse(p_break <0.001, "<0.001", sprintf("%.3f", p_break))
  
  newdat <- data.frame(log_gdppc = (1:70)*0.01+4.5,
                       pop50 = mean(sums$pop50))
  pred_prev <- predict.segmented(seg_fi, 
                                 newdata = newdat, 
                                 se.fit = TRUE)
  if (p_break>0.05/6){
    pred_prev <- predict(glm_fi, 
                         newdata = newdat, 
                         se.fit = TRUE)
  }
  pred_df <- data.frame(
    log_gdppc = (1:70)*0.01+4.5,
    `人均GDP` = 10^(log_gdppc),
    pred_fit = exp(pred_prev$fit)
  )
  
  plot_op <- ggplot() +
    geom_point(data = sums, aes(x = `人均GDP`, y = outcome, size = wgts), alpha = 0.5) +
    geom_line(aes(y = pred_df$pred_fit, x = pred_df$`人均GDP`, linetype = "dashed")) + 
    scale_x_log10() + 
    theme_bw() + 
    labs(
      x = "GDP per capita",
      y = "Prevalence") + 
    # no legend
    theme(legend.position = "none")
  return(list(plot_op, break_pt, p_op, sums, pred_df))
}

plot_hi <- plot_break(outcome_name = "hi_prev")
plot_vi <- plot_break(outcome_name = "vi_prev")
plot_sd <- plot_break(outcome_name = "sd_prev")
plot_ci <- plot_break(outcome_name = "ci_prev")
plot_dp <- plot_break(outcome_name = "dp_prev")
plot_ad <- plot_break(outcome_name = "ad_prev")


plot_all <- 
  (plot_vi[[1]] + ggtitle("A. Visual impairment")) + 
  (plot_hi[[1]] + ggtitle("B. Hearing impairment")) + 
  (plot_ci[[1]] + ggtitle("C. Cognitive impairment")) + 
  (plot_sd[[1]] + ggtitle("D. Sleep disorder")) + 
  (plot_dp[[1]] + ggtitle("E. Depressive symptoms")) + 
  (plot_ad[[1]] + ggtitle("F. ADL disability")) + 
  plot_layout(ncol = 2)
plot_all
ggsave("output/seg-prov.pdf", plot_all, width = 8, height = 12)

save.image("seg-prov.rda")


