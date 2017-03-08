
library(ggplot2)
library(BlandAltmanLeh)

df <- read.csv("Paired data large cluster.csv", sep = ';', header = T)

davidek <- as.numeric(unlist(df[1]))
pilot <- as.numeric(unlist(df[3]))

#Get data and stats:
ba.stats <- bland.altman.stats(davidek, pilot) #Always put in (Standard method, Method to compare)
means <- ba.stats$means
diffs <- ba.stats$diff
ba.data <- data.frame(means, diffs)

#Plot:
ggplot(ba.data,aes(means, diffs)) + 
  geom_point(alpha=0.15, color="darkcyan") + 
  coord_cartesian(xlim= c(60,220), ylim=c(-200,100)) + 
  geom_hline(yintercept = ba.stats$lines, col=c("darkgrey","black","darkgrey"), 
             linetype = 2, size=0.5) +
  scale_x_continuous(name="Mean of ratios") +
  scale_y_continuous(name="Difference in ratio")

#Save:
ggsave("BlandAltman.png", width = 10, height = 5)

#Statistics:
limits <- c(ba.stats$lower.limit, ba.stats$upper.limit)
inrangeDiffs <- subset(data.frame(diff=diffs), diff>limits[1] & diff<limits[2])
inrangePercentage <- 100*length(inrangeDiffs$diff)/length(diffs)
#Percentage data points within 1.96 std.dev:
print(inrangePercentage)
#Mean of differences, e.g. systematic error of method compared to standard:
print(ba.stats$mean.diffs)

