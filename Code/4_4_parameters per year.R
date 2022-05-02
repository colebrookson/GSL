pdf("Parameter distribution across the years - Atlantic cod.pdf")

# Fixed parameters
fixedMean <- sapply(model, function(x) x$summary.fixed[,1])
fixed.025 <- sapply(model, function(x) x$summary.fixed[,3])
fixed.975 <- sapply(model, function(x) x$summary.fixed[,5])

NameFixPara <- c("Intercept",
                 "Seal",
                 "Depth",
                 "Difference between surface\nand bottom temperature")

for(i in 1:nrow(fixedMean)){
  plot(1971:2020,
       fixedMean[i,],
       type = "l",
       xlab= "",
       ylab = "",
       main = NameFixPara[i],
       ylim = range(fixed.025[i,], fixed.975[i,]),
       lwd = 2)
  lines(1971:2020,fixed.025[i,], col = "grey")
  lines(1971:2020,fixed.975[i,], col = "grey")
  abline(v = 2018, col = "blue")
  abline(h = 0, col = "red")
}

# Hyperparameters
hyperMean <- sapply(model, function(x) x$summary.hyperpar[,1])
hyper.025 <- sapply(model, function(x) x$summary.hyperpar[,3])
hyper.975 <- sapply(model, function(x) x$summary.hyperpar[,5])

NameHyperPara <- c("Negative binomial standard deviation",
                 "Marginal standard deviation of the barrier model",
                 "Range of the barrier model")

for(i in 1:nrow(hyperMean)){
  plot(1971:2020,
       hyperMean[i,],
       type = "l",
       xlab= "",
       ylab = "",
       main = NameHyperPara[i],
       ylim = range(hyper.025[i,], hyper.975[i,]),
       lwd = 2,log = ifelse(i == 1, "y", ""))
  lines(1971:2020,hyper.025[i,], col = "grey")
  lines(1971:2020,hyper.975[i,], col = "grey")
  abline(v = 2018, col = "blue")
  abline(h = 0, col = "red")
}
 
dev.off()