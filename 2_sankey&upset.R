setwd('D:/Dataset/CHARLS')
rm(list = ls())

library(tidyverse)
library(tableone)
library(tableone)
library(ggalluvial)
library(ggsci)
library(ComplexUpset)


load('data/alldata_arfi_imputation_Vad.rda')

##1. Sankey plot ####
skdata <- alldata %>% 
  mutate(r1num = r1vi+r1hi+r1ci+r1sd+r1dp+r1ad, 
         r2num = r2vi+r2hi+r2ci+r2sd+r2dp+r2ad, 
         r3num = r3vi+r3hi+r3ci+r3sd+r3dp+r3ad, 
         r4num = r4vi+r4hi+r4ci+r4sd+r4dp+r4ad, 
         r5num = r5vi+r5hi+r5ci+r5sd+r5dp+r5ad, 
         r2num = ifelse(mo_oc==1&mo_t2<=2013, 'Death', r2num), 
         r3num = ifelse(mo_oc==1&mo_t2<=2015, 'Death', r3num), 
         r4num = ifelse(mo_oc==1&mo_t2<=2018, 'Death', r4num), 
         r5num = ifelse(mo_oc==1&mo_t2<=2020, 'Death', r5num), 
         r1num = factor(r1num), 
         r2num = factor(r2num), 
         r3num = factor(r3num), 
         r4num = factor(r4num), 
         r5num = factor(r5num)) %>% 
  filter(!is.na(r1num)) %>% 
  pivot_longer(cols = ends_with('num'),
               names_to = 'wave', 
               values_to = 'Status') %>% 
  mutate(
    wave = case_when(
      wave == 'r1num' ~ '2011', 
      wave == 'r2num' ~ '2013', 
      wave == 'r3num' ~ '2015', 
      wave == 'r4num' ~ '2018', 
      wave == 'r5num' ~ '2020', 
    ), 
    Status = case_when(
      Status == 0 ~ 'No functional impairments existing', 
      Status == 1 ~ 'One functional impairment existing', 
      Status == 2 ~ 'Two functional impairments coexisting', 
      Status == 3 ~ 'Three functional impairments coexisting', 
      Status == 4 ~ 'Four functional impairments coexisting', 
      Status == 5 ~ 'Five functional impairments coexisting', 
      Status == 6 ~ 'Six functional impairments coexisting', 
      TRUE ~ Status
    ), 
    wave = factor(wave), 
    Time = wave, 
    Status = factor(ifelse(is.na(Status), 'Not visited', Status), 
                    levels = c('No functional impairments existing', 
                               'One functional impairment existing', 
                               'Two functional impairments coexisting', 
                               'Three functional impairments coexisting', 
                               'Four functional impairments coexisting', 
                               'Five functional impairments coexisting', 
                               'Six functional impairments coexisting', 
                               'Not visited', 
                               'Death'), 
                    ordered = T)
  ) %>% 
  select(ID, wave, Status, Time)

summary(skdata)


skplot <- ggplot(skdata, aes(x = Time, stratum = Status, alluvium = ID, fill = Status, label = Status)) +
  geom_flow(color ="darkgray") + 
  geom_stratum() + labs(x = '', y = 'Number of participants', title = '') + 
  theme(legend.title = element_text(size=10), legend.text=element_text(size=8), legend.key.size = unit(35, "pt")) + 
  guides(shape = guide_legend(override.aes = list(size = 0.2))) + theme_minimal() +
  ggsci::scale_fill_lancet()

skplot

ggsave(file = 'output/skplot_Vad.pdf', plot = skplot, width = 11, height = 6)




skdata <- mutate(skdata, Status = ifelse(Status=='Death'|Status=='Not visited', NA, Status)-1)

table <- CreateTableOne(vars = 'Status', strata = 'Time', factorVars = 'Status', data = skdata)
tableout <- print(table, showAllLevels = T, formatOptions = list(big.mark = ','), quote = F, noSpaces = T, printToggle = F)
tableout
write.csv(tableout, file = 'output/sktable_Vad.csv')



##2. Upset plot ####
updata <- alldata %>% 
  mutate(r1num = r1vi+r1hi+r1ci+r1sd+r1dp+r1ad, 
         r2num = r2vi+r2hi+r2ci+r2sd+r2dp+r2ad, 
         r3num = r3vi+r3hi+r3ci+r3sd+r3dp+r3ad, 
         r4num = r4vi+r4hi+r4ci+r4sd+r4dp+r4ad, 
         r5num = r5vi+r5hi+r5ci+r5sd+r5dp+r5ad, 
         r2num = ifelse(mo_oc==1&mo_t2<=2013, 'Death', r2num), 
         r3num = ifelse(mo_oc==1&mo_t2<=2015, 'Death', r3num), 
         r4num = ifelse(mo_oc==1&mo_t2<=2018, 'Death', r4num), 
         r5num = ifelse(mo_oc==1&mo_t2<=2020, 'Death', r5num), 
         r1num = factor(r1num), 
         r2num = factor(r2num), 
         r3num = factor(r3num), 
         r4num = factor(r4num), 
         r5num = factor(r5num)) %>% 
  filter(!is.na(r1num))

summary(updata)

wave1 <- updata %>% 
  mutate(`Visual impairment`=r1vi, 
         `Hearing impairment`=r1hi, 
         `Cognitive impairment`=r1ci, 
         `Sleep disorder`=r1sd, 
         `Depressive symptoms`=r1dp, 
         `ADL disability`=r1ad) %>% 
  select(`Visual impairment`:`ADL disability`)


p1 <- upset(wave1, colnames(wave1), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=30,
                          themes=upset_modify_themes(
                            list('intersections_matrix'=theme(text=element_text(size=20)),
                                 "Intersection size"=theme(text=element_text(size=20)),
                                 'overall_sizes'=theme(text=element_text(size=20), axis.text.x=element_text(angle=90)))))


wave2 <- updata %>% 
  mutate(`Visual impairment`=r2vi, 
         `Hearing impairment`=r2hi, 
         `Cognitive impairment`=r2ci, 
         `Sleep disorder`=r2sd, 
         `Depressive symptoms`=r2dp, 
         `ADL disability`=r2ad) %>% 
  select(`Visual impairment`:`ADL disability`)


p2 <- upset(wave2, colnames(wave2), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=30,
            themes=upset_modify_themes(
              list('intersections_matrix'=theme(text=element_text(size=20)),
                   "Intersection size"=theme(text=element_text(size=20)),
                   'overall_sizes'=theme(text=element_text(size=20), axis.text.x=element_text(angle=90)))))


wave3 <- updata %>% 
  mutate(`Visual impairment`=r3vi, 
         `Hearing impairment`=r3hi, 
         `Cognitive impairment`=r3ci, 
         `Sleep disorder`=r3sd, 
         `Depressive symptoms`=r3dp, 
         `ADL disability`=r3ad) %>% 
  select(`Visual impairment`:`ADL disability`)


p3 <- upset(wave3, colnames(wave3), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=30,
            themes=upset_modify_themes(
              list('intersections_matrix'=theme(text=element_text(size=20)),
                   "Intersection size"=theme(text=element_text(size=20)),
                   'overall_sizes'=theme(text=element_text(size=20), axis.text.x=element_text(angle=90)))))


wave4 <- updata %>% 
  mutate(`Visual impairment`=r4vi, 
         `Hearing impairment`=r4hi, 
         `Cognitive impairment`=r4ci, 
         `Sleep disorder`=r4sd, 
         `Depressive symptoms`=r4dp, 
         `ADL disability`=r4ad) %>% 
  select(`Visual impairment`:`ADL disability`)


p4 <- upset(wave4, colnames(wave4), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=30,
            themes=upset_modify_themes(
              list('intersections_matrix'=theme(text=element_text(size=20)),
                   "Intersection size"=theme(text=element_text(size=20)),
                   'overall_sizes'=theme(text=element_text(size=20), axis.text.x=element_text(angle=90)))))

wave5 <- updata %>% 
  mutate(`Visual impairment`=r5vi, 
         `Hearing impairment`=r5hi, 
         `Cognitive impairment`=r5ci, 
         `Sleep disorder`=r5sd, 
         `Depressive symptoms`=r5dp, 
         `ADL disability`=r5ad) %>% 
  select(`Visual impairment`:`ADL disability`)


p5 <- upset(wave5, colnames(wave5), name='Co-occurrence of ARFIs', width_ratio=0.25, min_size=30,
            themes=upset_modify_themes(
              list('intersections_matrix'=theme(text=element_text(size=20)),
                   "Intersection size"=theme(text=element_text(size=20)),
                   'overall_sizes'=theme(text=element_text(size=20), axis.text.x=element_text(angle=90)))))


ggsave(file = 'output/upset_p1_Vad.pdf', plot = p1, width = 16, height = 9.6)
ggsave(file = 'output/upset_p2_Vad.pdf', plot = p2, width = 16, height = 9.6)
ggsave(file = 'output/upset_p3_Vad.pdf', plot = p3, width = 16, height = 9.6)
ggsave(file = 'output/upset_p4_Vad.pdf', plot = p4, width = 16, height = 9.6)
ggsave(file = 'output/upset_p5_Vad.pdf', plot = p5, width = 16, height = 9.6)


wave_data <- wave5

pattern5 <- as.data.frame(table(
  wave_data$`Visual impairment`, 
  wave_data$`Hearing impairment`, 
  wave_data$`Cognitive impairment`, 
  wave_data$`Sleep disorder`, 
  wave_data$`Depressive symptoms`, 
  wave_data$`ADL disability`
))

write.csv(pattern, 'output/copattern1.csv', row.names = F)
write.csv(pattern, 'output/copattern2.csv', row.names = F)
write.csv(pattern, 'output/copattern3.csv', row.names = F)
write.csv(pattern, 'output/copattern4.csv', row.names = F)
write.csv(pattern, 'output/copattern5.csv', row.names = F)

allpattern <- left_join(pattern1, pattern2, by = names(pattern1)[1:6]) %>% 
  left_join(pattern3, by = names(pattern1)[1:6]) %>% 
  left_join(pattern4, by = names(pattern1)[1:6]) %>% 
  left_join(pattern5, by = names(pattern1)[1:6])

write.csv(allpattern, 'output/allpattern.csv', row.names = F)

