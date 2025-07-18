setwd('~')
rm(list = ls())

library(tidyverse)
library(tableone)
library(survival)
library(pheatmap)
library(forestploter)
library(grid)


# dyads-mortality and heatmap ####
### dyads
arfis <- c('r1vi', 'r1hi', 'r1ci', 'r1sd', 'r1dp', 'r1ad')
pair_HR <- matrix(data = NA, nrow = 5, ncol = 5)

for (i in 1:5) {

  row_HR <- rep('', 5)
  
  for (j in (i+1):6){
    
    pairdata <- alldata %>% 
      mutate(e1 = alldata[[arfis[i]]], 
             e2 = alldata[[arfis[j]]], 
             exposure = e1+e2, 
             OtherARFIs = rowSums(alldata[arfis[-c(i,j)]], na.rm = T), 
             OtherARFIs = ifelse(OtherARFIs>0, 1, OtherARFIs))
    
    pair_cox <- coxph(Surv(mo_t2-2011, mo_oc)~factor(exposure)+
                        AGE+GENDER+EDU+MARR+HUKOU+SET+INCOME+BMI+SMK+DRK+STR+CAN+HBP+DIA+HEA+OtherARFIs, data=pairdata)
    
    tmpdata <- filter(pairdata, !is.na(exposure), !is.na(mo_oc))
    
    message(sum(tmpdata$mo_oc), '/', sum(tmpdata$mo_t2-2011))
    
    row_HR[j-1] <- exp(pair_cox$coefficients[2])
    print(arfis[i])
    print(arfis[j])
    print(ShowRegTable(pair_cox, printToggle = F)[2, ])
    
  }
  
  pair_HR[i, ] <- row_HR
  
}


pair_HR <- t(pair_HR)
colnames(pair_HR) <- c('Visual impairment', 'Hearing impairment', 'Cognitive impairment', 'Sleep disorder', 'Depressive symptons')
rownames(pair_HR) <- c('Hearing impairment', 'Cognitive impairment', 'Sleep disorder', 'Depressive symptons', 'ADL disability')

pair_HR

write.csv(x = pair_HR, file = 'output/dyads_HR_Vad.csv')


### heatmap
hr_m <- read.csv('output/dyads_HR_Vad.csv', row.names = 1)
hr_m[1, 2] <- 0
heatplot <- pheatmap(hr_m, cluster_row = F, cluster_cols = F, angle_col = 315,
                     display_numbers = T, 
                     border='white', color = colorRampPalette(c('#6895c7', "white","#684f78"))(100))

ggsave(filename = 'output/Hheatplot_Vad.pdf', plot = heatplot, height = 3.5, width = 5)


numberdata <- alldata %>% 
  mutate(number = rowSums(alldata[c('r1vi', 'r1hi', 'r1ci', 'r1sd', 'r1dp', 'r1ad')], na.rm = T), 
         number = ifelse(is.na(r1vi)&is.na(r1hi)&is.na(r1ci)&is.na(r1sd)&is.na(r1dp)&is.na(r1pf), NA, number), 
         number_c = factor(ifelse(number>4, 4, number), 
                           labels = c('Without any ARFI', 
                                      'With one ARFI', 
                                      'With two ARFIs', 
                                      'With three ARFIs', 
                                      'With more than three ARFIs')), 
         mo_t2 = mo_t2-2011)
summary(numberdata)


age_number <- numberdata %>% 
  mutate(number = ifelse(number>=4, 4, number)) %>% 
  filter(!is.na(number)) %>% 
  group_by(AGC) %>% 
  summarize(n = n())

age_number_specific <- numberdata %>% 
  mutate(number = ifelse(number>=4, 4, number)) %>% 
  filter(!is.na(number)) %>% 
  group_by(number, AGC) %>% 
  summarize(n = n())

age_number_specific$pre <- age_number_specific$n/age_number$n
age_number_specific$totalnumber <- rep(age_number$n, 5)
age_number_specific$lower <- age_number_specific$pre - 1.96 * sqrt(age_number_specific$pre*(1-age_number_specific$pre)/age_number_specific$totalnumber)
age_number_specific$upper <- age_number_specific$pre + 1.96 * sqrt(age_number_specific$pre*(1-age_number_specific$pre)/age_number_specific$totalnumber)

write.csv(age_number_specific, file = 'output/age_number_specific.csv')

num_cox <- coxph(Surv(mo_t2-2011, mo_oc)~number_c+
                   AGE+GENDER+EDU+INCOME+BMI+SMK+DRK+STR+CAN+HBP+DIA+HEA, data=numberdata)
each_cox <- coxph(Surv(mo_t2-2011, mo_oc)~number+
                    AGE+GENDER+EDU+INCOME+BMI+SMK+DRK+STR+CAN+HBP+DIA+HEA, data=numberdata)
summary(num_cox)
summary(each_cox)


forest_results <- data.frame(rbind(c(NA,NA,NA),
                                   c(1,1,1),
                                   cbind(summary(num_cox)[["conf.int"]][1:4, c(1,3,4)]),
                                   cbind(summary(each_cox)[["conf.int"]][, c(1,3,4)]))[1:7,])
colnames(forest_results) <- c('mean', 'lower', 'upper')
rownames(forest_results) <- NULL

fordata <- data.frame(
  numberdata %>%
    filter(!is.na(number_c)) %>% 
    mutate(t = mo_t2) %>%
    group_by(number_c) %>%
    summarize(sum(mo_oc, na.rm = T), sum(t, na.rm = T))
)
fordata

c(ShowRegTable(num_cox, printToggle = F)[1:4, 2], ShowRegTable(each_cox, printToggle = F)[1, 2])


tabletext <- cbind(
  c('','Without any ARFI', 'With one ARFI', 'With two ARFIs', 'With three ARFIs', 'With more than three ARFIs', 'Each additional ARFI'),
  paste0(c('Cases', fordata[,2], 7981),'/',c('Person-years', fordata[,3], 227657)),
  c('HR [95% CI]',"1.00 [Reference]",
    as.vector(ShowRegTable(num_cox, printToggle = F)[1,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[2,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[3,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[4,1]),
    as.vector(ShowRegTable(each_cox, printToggle = F)[1,1])),
  c(paste('P-value'), '', c(ShowRegTable(num_cox, printToggle = F)[1:4, 2], ShowRegTable(each_cox, printToggle = F)[1, 2]))
)
rownames(tabletext) <- NULL


specific_results <- as.data.frame(
  cbind(
    c('Specific ARFI', 'Visual impairment', 'Hearing impairment', 'Cognitive impairment', 
      'Sleep disorder', 'Depressive symptoms', 'ADL disability', 'The number of ARFIs'), 
    c(
      ' ', 
      paste0(nrow(filter(alldata, r1vi==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1vi==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1hi==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1hi==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1ci==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1ci==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1sd==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1sd==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1dp==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1dp==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1ad==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1ad==0, mo_t2>2011)$mo_t2-2011)), 
      ' '
    ), 
    c(' ', 
      '1.07 [0.93, 1.23]', '0.95 [0.85, 1.07]', '1.41 [1.21, 1.64]', 
      '0.96 [0.86, 1.07]', '1.20 [1.08, 1.33]', '1.64 [1.48, 1.80]', 
      ' '), 
    rep(NA, 8), 
    c(NA, 1.07, 0.95, 1.41, 0.96, 1.20, 1.64, NA), 
    c(NA, 0.93, 0.85, 1.21, 0.86, 1.08, 1.48, NA), 
    c(NA, 1.23, 1.07, 1.64, 1.07, 1.33, 1.80, NA)
  )
)

df1 <- cbind(tabletext, forest_results)
names(df1) <- c(df1[1,1:4], 'mean', 'lower', 'upper')
names(specific_results) <- c(df1[1,1:4], 'mean', 'lower', 'upper')
df1 <- df1[-1,]
df1 <- rbind(specific_results, df1)
df1$mean <- as.numeric(df1$mean)
df1$lower <- as.numeric(df1$lower)
df1$upper <- as.numeric(df1$upper)
df1$` ` <- paste(rep(" ", 20), collapse = " ")


tm <- forest_theme(
  base_size = 10,
  ci_pch = 16,
  ci_col = '#6d4e84',
  ci_lty = 1, 
  ci_lwd = 1.8,
  ci_Theight = 0.3,
  refline_lwd = 1.5,
  refline_col = "grey60",
  core=list(bg_params=list(fill =c("#fefefe", "#ede7ef")))
)

p1 <- forest(df1[2:7, c(1, 2, 8, 3)],
             est = df1[2:7, ]$mean,
             lower = df1[2:7, ]$lower,
             upper = df1[2:7, ]$upper,
             ci_column = 3, 
             xlim = c(0.8, 2),
             ticks_at = c(1, 1.5, 2),
             ref_line = 1,
             theme = tm,
             title = '') %>% 
  add_border(., part = "header", where = "bottom") %>% 
  add_border(., part = "header", where = "top") 
p1

ggsave(filename = 'output/figure3a.pdf', plot = p1, width = 10, height = 5)


p2 <- forest(df1[9:14, c(1, 2, 8, 3)],
             est = df1[9:14, ]$mean,
             lower = df1[9:14, ]$lower,
             upper = df1[9:14, ]$upper,
             ci_column = 3, 
             xlim = c(0.8, 2),
             ticks_at = c(1, 1.5, 2),
             ref_line = 1,
             theme = tm,
             title = '') %>% 
  add_border(., part = "header", where = "bottom") %>% 
  add_border(., part = "header", where = "top")

p2

ggsave(filename = 'output/figure4a.pdf', plot = p2, width = 10, height = 5)









cor_data <- alldata %>% 
  select(ID, AGE, GENDER, EDU, INCOME, PA, BMI, SMK, DRK, STR, CAN, MEM, PSY, HBP, DIA, HEA)
load('data/alldata_arfi_no_permanet.rda')

# alldata <- alldata %>% 
#   mutate(changes_c = as.numeric(change1115)-5, 
#          changes = ifelse(changes_c<=0, 0, ifelse(changes_c>3, 3, changes_c)), 
#          changes = factor(changes)) %>% 
#   select(ID, mo_t2, mo_oc, changes, changes_c) %>% 
#   left_join(cor_data, by = 'ID')
# 
# change_cox <- coxph(Surv(mo_t2-2011, mo_oc)~changes+
#                       AGE+GENDER+EDU+INCOME+PA+BMI+SMK+DRK+STR+CAN+MEM+PSY+HBP+DIA+HEA, data=alldata)
# each_cox <- coxph(Surv(mo_t2-2011, mo_oc)~changes_c+
#                     AGE+GENDER+EDU+INCOME+PA+BMI+SMK+DRK+STR+CAN+MEM+PSY+HBP+DIA+HEA, data=alldata)


dt <- alldata %>% 
  mutate(mo_t2 = mo_t2-2015, 
         change1115 = as.numeric(change1115)-5, 
         r1vi = 0, 
         r1hi = 0, 
         baselinenumber = r1vi+r1hi+r1ci+r1sd+r1dp+r1ad,
         change1115_2 = ifelse(change1115<0, 0, change1115),
         change1115_c = factor(ifelse(change1115>3, 3, ifelse(change1115<0, 0, change1115)))) %>% 
  filter(mo_t2>0) %>% 
  select(ID, mo_t2, mo_oc, change1115, change1115_2, change1115_c, baselinenumber) %>% 
  left_join(cor_data, by = 'ID')

summary(dt)

dd <- datadist(dt) 
options(datadist='dd')
S <- Surv(dt$mo_t2,dt$mo_oc==1)
fit <- cph(S ~ rcs(change1115, 4)+AGE+GENDER+EDU+INCOME+PA+BMI+SMK+DRK+STR+CAN+MEM+PSY+baselinenumber, x=TRUE, y=TRUE,data=dt)

Pre0 <-rms::Predict(fit, change1115, fun=exp, type="predictions", ref.zero=TRUE, conf.int = 0.95, digits=2)

Pre0 <- Pre0 %>% 
  filter(change1115>=-2, change1115<=3)

rcsplot <- ggplot()+
  geom_line(data=Pre0, aes(change1115, yhat), alpha=0.7)+
  geom_ribbon(data=Pre0, aes(change1115, ymin=lower, ymax=upper), alpha=0.1)+
  geom_hline(yintercept=1, linetype=2, size=0.75) +
  # coord_cartesian(ylim = c(0, 3.5)) + 
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  labs(x = 'Changes in number of ARFIs from 2011 to 2015', y = 'HR (95% CI)', title = '') +
  theme_classic()
rcsplot
ggsave(plot = rcsplot, filename = 'output/rcsplot.pdf', width = 6, height = 4.2)


changemodel <- coxph(Surv(mo_t2, mo_oc)~
                       change1115_c+AGE+GENDER+EDU+INCOME+PA+BMI+SMK+DRK+STR+CAN+MEM+PSY+HBP+DIA+HEA+baselinenumber, data = dt)
ShowRegTable(changemodel)

changemodel2 <- coxph(Surv(mo_t2, mo_oc)~
                        change1115_2+AGE+GENDER+EDU+INCOME+PA+BMI+SMK+DRK+STR+CAN+MEM+PSY+HBP+DIA+HEA+baselinenumber, data = dt)
ShowRegTable(changemodel2)


paste0(nrow(filter(dt, change1115_c==0, mo_oc==1, mo_t2>0)), '/', sum(filter(dt, change1115_c==0, mo_t2>0)$mo_t2))
paste0(nrow(filter(dt, change1115_c==1, mo_oc==1, mo_t2>0)), '/', sum(filter(dt, change1115_c==1, mo_t2>0)$mo_t2))
paste0(nrow(filter(dt, change1115_c==2, mo_oc==1, mo_t2>0)), '/', sum(filter(dt, change1115_c==2, mo_t2>0)$mo_t2))
paste0(nrow(filter(dt, change1115_c==3, mo_oc==1, mo_t2>0)), '/', sum(filter(dt, change1115_c==3, mo_t2>0)$mo_t2))
paste0(nrow(filter(dt, mo_oc==1, mo_t2>0)), '/', sum(filter(dt, mo_t2>0)$mo_t2))



df1[9:13, 1] <- c('no increment', '1 increment', '2 increment', '>=3 increment', 'each additional increment')
df1[9:13, 2] <- c("201/14025", "67/4031", "32/1515", "17/600", "317/20171")
df1[10:13, 3] <- c('1.12 [0.85, 1.49]', '1.54 [1.05, 2.26]', '1.91 [1.15, 3.17]', '1.21 [1.07, 1.37]')
df1[10:13, 5] <- c(1.12, 1.54, 1.91, 1.21)
df1[10:13, 6] <- c(0.85, 1.05, 1.15, 1.07)
df1[10:13, 7] <- c(1.49, 2.26, 3.17, 1.37)

p3 <- forest(df1[9:13, c(1, 2, 8, 3)],
             est = df1[9:13, ]$mean,
             lower = df1[9:13, ]$lower,
             upper = df1[9:13, ]$upper,
             ci_column = 3, 
             # xlim = c(0.8, 2),
             # ticks_at = c(1, 1.5, 2),
             ref_line = 1,
             theme = tm,
             title = '') %>% 
  add_border(., part = "header", where = "bottom") %>% 
  add_border(., part = "header", where = "top") # %>% 
# edit_plot(row = c(1, 8), gp = gpar(fontface = "bold"))
p3

ggsave(filename = 'output/figure4b.pdf', plot = p3, width = 10, height = 5)



# number of ARFI-mortality and forest (Three version)####
numberdata <- alldata %>% 
  mutate(number = rowSums(alldata[c('r1ci', 'r1dp', 'r1ad')], na.rm = T), 
         number = ifelse(is.na(r1ci)&is.na(r1dp)&is.na(r1pf), NA, number), 
         number_c = factor(number,
                           labels = c('Without any ARFI', 
                                      'With one ARFI', 
                                      'With two ARFIs', 
                                      'With three ARFIs')), 
         mo_t2 = mo_t2-2011)
# summary(numberdata)


age_number <- numberdata %>% 
  filter(!is.na(number)) %>% 
  group_by(AGC) %>% 
  summarize(n = n())

age_number_specific <- numberdata %>% 
  filter(!is.na(number)) %>% 
  group_by(number, AGC) %>% 
  summarize(n = n())

age_number_specific$pre <- age_number_specific$n/age_number$n
age_number_specific$totalnumber <- rep(age_number$n, 4)
age_number_specific$lower <- age_number_specific$pre - 1.96 * sqrt(age_number_specific$pre*(1-age_number_specific$pre)/age_number_specific$totalnumber)
age_number_specific$upper <- age_number_specific$pre + 1.96 * sqrt(age_number_specific$pre*(1-age_number_specific$pre)/age_number_specific$totalnumber)

write.csv(age_number_specific, file = 'output/age_number_specific.csv')

num_cox <- coxph(Surv(mo_t2-2011, mo_oc)~number_c+
                   AGE+GENDER+EDU+MARR+HUKOU+SET+INCOME+BMI+SMK+DRK+STR+CAN+HBP+DIA+HEA, data=numberdata)
each_cox <- coxph(Surv(mo_t2-2011, mo_oc)~number+
                    AGE+GENDER+EDU+MARR+HUKOU+SET+INCOME+BMI+SMK+DRK+STR+CAN+HBP+DIA+HEA, data=numberdata)
summary(num_cox)
summary(each_cox)


forest_results <- data.frame(rbind(c(NA,NA,NA),
                                   c(1,1,1),
                                   cbind(summary(num_cox)[["conf.int"]][1:3, c(1,3,4)]),
                                   cbind(summary(each_cox)[["conf.int"]][, c(1,3,4)]))[1:6,])
colnames(forest_results) <- c('mean', 'lower', 'upper')
rownames(forest_results) <- NULL

fordata <- data.frame(
  numberdata %>%
    filter(!is.na(number_c)) %>% 
    mutate(t = mo_t2) %>%
    group_by(number_c) %>%
    summarize(sum(mo_oc, na.rm = T), sum(t, na.rm = T))
)
fordata

c(ShowRegTable(num_cox, printToggle = F)[1:4, 2], ShowRegTable(each_cox, printToggle = F)[1, 2])


tabletext <- cbind(
  c('','Without any ARFI', 'With one ARFI', 'With two ARFIs', 'With three ARFIs', 'Each additional ARFI'),
  paste0(c('Cases', fordata[,2], 7981),'/',c('PY', fordata[,3], 227657)),
  c('HR [95% CI]',"1.00 [Reference]",
    as.vector(ShowRegTable(num_cox, printToggle = F)[1,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[2,1]),
    as.vector(ShowRegTable(num_cox, printToggle = F)[3,1]),
    # as.vector(ShowRegTable(num_cox, printToggle = F)[4,1]),
    # as.vector(ShowRegTable(num_cox, printToggle = F)[5,1]),
    as.vector(ShowRegTable(each_cox, printToggle = F)[1,1])),
  c(paste('P-value'), '', c(ShowRegTable(num_cox, printToggle = F)[1:3, 2], ShowRegTable(each_cox, printToggle = F)[1, 2]))
)
rownames(tabletext) <- NULL

# pdf('forest.pdf', width = 10, height = 3.5)
# forestplot(labeltext = tabletext, forest_results, graph.pos=4, graphwidth = unit(7, "cm"), lineheight = unit(0.8, 'cm'),
#            hrzl_lines = list('1' = gpar(lwd=2),
#                              '2' = gpar(lwd=1)), new_page = T, align = 'l', colgap = unit(0.01, "npc"), boxsize = 0.2, 
#            is.summary=c(T, rep(F,7)), xlog=T, vertices = T, title = 'The association of the number of ARFIs with mortality')
# dev.off()

specific_results <- as.data.frame(
  cbind(
    c('Specific ARFI', 'Visual impairment', 'Hearing impairment', 'Cognitive impairment', 
      'Sleep disorder', 'Depressive symptoms', 'ADL disability', 'The number of key ARFIs'), 
    c(
      ' ', 
      paste0(nrow(filter(alldata, r1vi==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1vi==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1hi==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1hi==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1ci==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1ci==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1sd==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1sd==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1dp==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1dp==0, mo_t2>2011)$mo_t2-2011)), 
      paste0(nrow(filter(alldata, r1ad==0, mo_oc==1, mo_t2>2011)), '/', sum(filter(alldata, r1ad==0, mo_t2>2011)$mo_t2-2011)), 
      ' '
    ), 
    c(' ', 
      '1.08 [0.94, 1.24]', '0.96 [0.85, 1.07]', '1.40 [1.20, 1.63]', 
      '0.97 [0.87, 1.08]', '1.18 [1.06, 1.30]', '1.66 [1.51, 1.83]', 
      ' '), 
    rep(NA, 8), 
    c(NA, 1.08, 0.96, 1.40, 0.97, 1.18, 1.66, NA), 
    c(NA, 0.94, 0.85, 1.20, 0.87, 1.06, 1.51, NA), 
    c(NA, 1.24, 1.07, 1.63, 1.08, 1.30, 1.83, NA)
  )
)

df1 <- cbind(tabletext, forest_results)
names(df1) <- c(df1[1,1:4], 'mean', 'lower', 'upper')
names(specific_results) <- c(df1[1,1:4], 'mean', 'lower', 'upper')
df1 <- df1[-1,]
df1 <- rbind(specific_results, df1)
df1$mean <- as.numeric(df1$mean)
df1$lower <- as.numeric(df1$lower)
df1$upper <- as.numeric(df1$upper)
df1$` ` <- paste(rep(" ", 20), collapse = " ")


tm <- forest_theme(
  base_size = 10,
  ci_pch = 16,
  ci_col = '#6d4e84',
  ci_lty = 1, 
  ci_lwd = 1.8,
  ci_Theight = 0.3,
  refline_lwd = 1.5,
  refline_col = "grey60",
  core=list(bg_params=list(fill =c("#fefefe", "#ede7ef")))
)

# p1 <- forest(df1[2:7, c(1, 2, 8, 3)],
#              est = df1[2:7, ]$mean,
#              lower = df1[2:7, ]$lower,
#              upper = df1[2:7, ]$upper,
#              ci_column = 3, 
#              xlim = c(0.8, 2),
#              ticks_at = c(1, 1.5, 2),
#              ref_line = 1,
#              theme = tm,
#              title = '') %>% 
#   add_border(., part = "header", where = "bottom") %>% 
#   add_border(., part = "header", where = "top") # %>% 
# # edit_plot(row = c(1, 8), gp = gpar(fontface = "bold"))
# p1

p1 <- forest(df1[, c(1, 2, 8, 3)],
             est = df1$mean,
             lower = df1$lower,
             upper = df1$upper,
             ci_column = 3, 
             xlim = c(0.8, 2),
             x_trans = 'log10', 
             ticks_at = c(1, 1.5, 2),
             ref_line = 1,
             theme = tm,
             title = '') %>% 
  add_border(., part = "header", where = "bottom") %>% 
  add_border(., part = "header", where = "top") %>% 
  edit_plot(row = c(1, 8), gp = gpar(fontface = "bold"))
p1


ggsave(filename = 'output/figure3a_0810.pdf', plot = p1, width = 7.5, height = 4.5)


p2 <- forest(df1[9:13, c(1, 2, 8, 3)],
             est = df1[9:13, ]$mean,
             lower = df1[9:13, ]$lower,
             upper = df1[9:13, ]$upper,
             ci_column = 3, 
             xlim = c(0.8, 2),
             ticks_at = c(1, 1.5, 2),
             ref_line = 1,
             theme = tm,
             title = '') %>% 
  add_border(., part = "header", where = "bottom") %>% 
  add_border(., part = "header", where = "top") # %>% 
# edit_plot(row = c(1, 8), gp = gpar(fontface = "bold"))
p2

ggsave(filename = 'output/figure4a_threeversion.pdf', plot = p2, width = 10, height = 5)





