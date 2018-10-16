dat <- read.csv("Dropbox/Tanja&XinW/Rfiles/CBD-O/FullObeseYear1.csv")
library(ggplot2)
pdf("Research/Obesity/docs/Ageplot.pdf",width = 6,height = 4)
ggplot(data = dat,aes(x = IYEAR, y = PropObese, group = AGE, color = AGE)) + 
  geom_point() + geom_line() + xlab("Year") + ylab("Proportion of obesity") +
  theme_bw()
dev.off()
