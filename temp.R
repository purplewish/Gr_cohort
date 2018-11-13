res_year1 <- refit_cohort(year = year, age = age, y = dat$PropObese, x = x0,group = c(1:27,27),model="year") 

res_gr_year1 <- refit_group(year = year, age = age, y = dat$PropObese, x = x0,group = c(1:27,27),model="year")

res_c <- refit_cohort2(year = year, age = age,y = dat$PropObese, x = x0,group.individual = groupest,group.cohort = groupestc,model = "year")

preddat1 <- dat
preddat1$est <- res_year1$estimates
preddat1$group <- as.factor(rep(c(1:27,27), each = length(unique(age))))
preddat1$curve <- res_year1$curve
preddat1$curve_gr <- res_gr_year1$estimates

ggplot(data = preddat1) + 
  geom_point(aes(x = AGE, y = PropObese, group = group, color = group)) + 
  #geom_line(aes(x = AGE, y = curve, group = group, color = group)) +
  geom_line(aes(x = AGE, y = curve_gr, group = group, color = group)) +
  theme_bw() 