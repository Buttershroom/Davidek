

df <- read.csv("Ttest.csv", sep = ';', header = T)

a <- as.numeric(unlist(df[1:3752, 1]))
b <- as.numeric(unlist(df[1:3752, 2]))

a = log(a)
b = log(b)
max(a)
max(b)
cmax <- max(max(a), max(b))
cmin <- min(min(a), min(b))

cmin = -3
cmax = 5

#Reset plot:
dev.off()

# Kernel Density Plot

da <- density(a) # returns the density data
plot(da, xlim = c(cmin, cmax), col="red", ylim = c(0, 0.6), lwd=2.5, xlab='') # plots the results
lines(density(b, na.rm = T), col = "blue", lwd=2.5)

legend('topright', c('Davidek', 'Pilot'), lty=c(1,1), lwd=c(2.5,2.5),col=c('red','blue'))

#Histograms:
hist(a, xlim = c(cmin, cmax), col = "red", breaks=400)
lines(density(a, na.rm = T), col = "red")
hist(b, xlim = c(cmin, cmax), add = T, col = "blue", breaks = 100)
lines(density(b, na.rm = T), col = "blue")
hist(c(0), xlim = c(cmin, cmax), add = F, col = "blue", ylim = c(0, 1))
median(a, na.rm = T)
median(b, na.rm = T)
t.test(a,b)

#CVs:
cv <- function(x){sd(x, na.rm = T)/mean(x, na.rm = T)}
cv(a)
cv(b)

#Boxplots:
boxplot(a,b)
boxplot(b)

#
