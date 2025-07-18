setwd('~')
rm(list = ls())

library(tidyverse)
library(tableone)
library(survey)
library(patchwork)
library(survival)


load('data/alldata_arfi_imputation_Vad.rda')

## Incidence rate ####
ir <- function(arfi, arfi_name, df = alldata) {
  
  d1 <- df[c('weight', 'AGC', paste0('r1', arfi), paste0('r2', arfi))]
  d2 <- df[c('weight', 'AGC', paste0('r2', arfi), paste0('r3', arfi))]
  names(d1) <- c('weight', 'AGC', 'arfi1', 'arfi2')
  names(d2) <- c('weight', 'AGC', 'arfi1', 'arfi2')
  newdata <- rbind(d1, d2)
  
  casedata <- df[c('AGC', paste0('r1', arfi), paste0('r2', arfi), paste0('r3', arfi))]
  names(casedata) <- c('AGC', 'arfi1', 'arfi2', 'arfi3')
  
  case_all <- casedata %>%
    filter(arfi1==0) %>%
    group_by(AGC) %>%
    summarise(case_all = paste0(sum(arfi3, na.rm = T), '/', sum(!is.na(arfi3))))
  
  ir_data <- newdata %>% 
    filter(arfi1==0, !is.na(arfi2)) %>% 
    group_by(AGC) %>%
    summarise(incidence_rate = sum(arfi2*weight)/(2*sum(weight))*1000) %>% 
    left_join(case_all, by = 'AGC') %>%
    mutate(ARFI = arfi_name) %>% 
    select(ARFI, AGC, case_all, incidence_rate)
  
  return(ir_data)
  
}


ir_all <- rbind(
  ir('vi', 'Visual impairment'),
  ir('hi', 'Hearing impairment'),
  ir('ci', 'Cognitive impairment'),
  ir('sd', 'Sleep disorder'),
  ir('dp', 'Depressive symptons'),
  ir('ad', 'ADL disability')
)

ir_male <- rbind(
  ir('vi', 'Visual impairment', filter(alldata, GENDER==1)),
  ir('hi', 'Hearing impairment', filter(alldata, GENDER==1)),
  ir('ci', 'Cognitive impairment', filter(alldata, GENDER==1)),
  ir('sd', 'Sleep disorder', filter(alldata, GENDER==1)),
  ir('dp', 'Depressive symptons', filter(alldata, GENDER==1)),
  ir('ad', 'ADL disability', filter(alldata, GENDER==1))
)

ir_female <- rbind(
  ir('vi', 'Visual impairment', filter(alldata, GENDER==2)),
  ir('hi', 'Hearing impairment', filter(alldata, GENDER==2)),
  ir('ci', 'Cognitive impairment', filter(alldata, GENDER==2)),
  ir('sd', 'Sleep disorder', filter(alldata, GENDER==2)),
  ir('dp', 'Depressive symptons', filter(alldata, GENDER==2)),
  ir('ad', 'ADL disability', filter(alldata, GENDER==2))
)

figure1 <- ggplot(data = ir_all, aes(x = AGC, y = incidence_rate, group = ARFI, color = ARFI, shape = ARFI)) + 
  geom_point() + geom_line(linewidth = 0.8) + theme_classic() + 
  theme(legend.position = 'none') + 
  labs(x = 'Age groups', y = 'Incidence rate (case/1000 person years)', title = 'A. Both genders') + 
  coord_cartesian(ylim = c(0,200)) +
  ggsci::scale_color_lancet()

figure1_male <- ggplot(data = ir_male, aes(x = AGC, y = incidence_rate, group = ARFI, color = ARFI, shape = ARFI)) + 
  geom_point() + geom_line(linewidth = 0.8) + theme_classic() + 
  theme(legend.title=element_blank()) + 
  labs(x = 'Age groups', y = 'Incidence rate (case/1000 person years)', title = 'C. Males') + 
  coord_cartesian(ylim = c(0,200)) + 
  ggsci::scale_color_lancet()

figure1_female <- ggplot(data = ir_female, aes(x = AGC, y = incidence_rate, group = ARFI, color = ARFI, shape = ARFI)) + 
  geom_point() + geom_line(linewidth = 0.8) + theme_classic() + 
  theme(legend.title=element_blank()) + 
  labs(x = 'Age groups', y = 'Incidence rate (case/1000 person years)', title = 'B. Females') + 
  coord_cartesian(ylim = c(0,200)) + 
  ggsci::scale_color_lancet()



p1 <- figure1 + (figure1_female / figure1_male) + plot_layout(guides = 'collect', widths = c(2, 1))
ggsave(file = 'output/figure1_Vad.pdf', plot = p1, width = 12, height = 6.4) #8.5 6.4
