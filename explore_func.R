library(infotheo)
library(arules)
library(nlme)
library(mgcv)
library(scam)

# required local variables: marr and rseqt
.plot_gene <- function(gene) {
  r <- cor(marrt[gene, ], rseqt[gene, ], use="complete.obs");
  plot(
    marrt[gene, ], rseqt[gene, ],
    xlab = "microarray log expression",
    ylab = "RNA-seq log expression",
    xlim = marrt.lim,
    ylim = rseq.lim,
    main = sprintf("%s (r = %s, r^2 = %s)", gene,
                   format(r, digits=2), format(r^2, digits=2))
  )
  r
}


.lm_gene <- function(gene) {
  lm(rseqt[gene, ] ~ marrt[gene, ])
}


# Automatically generate scatter plots according to pic_num and probe_list
scatter_plot <- function( pic_num, non_func_list) {
  stopifnot(pic_num > 0 & pic_num < dim(non_func_list)[1])
  for (i in 1:pic_num) {
    print(.plot_gene(non_func_list$gene[i]))
    readline("Press Enter to continue...")
  }
}



# Calculate variance explained of first component by PCA for each probes
pca.ratio <- function(mar ,rna ){
  pca_1 <- c()
  for (i in 1:nrow(mar)) {
    gene_i <- cbind(mar[i, ],rna[i, ])
    pca_i <- princomp(gene_i)
    sdev_first_pc <- pca_i$sdev[1]
    var_first <- (sdev_first_pc^2) / sum(pca_i$sdev^2)
    pca_1 <- append(pca_1, var_first)
  }
  return(pca_1)
}

# fraction of explained variance
fev.func<- function(x, x.hat) {
  ssr <- sum((x - x.hat)^2);
  sst <- sum((x - mean(x))^2);
  1 - ssr / (sst+1e-8)
}


# Calculate mutualinformation for each probes between Microarray and RNA-seq
mutual_info <- function(mar = marrt, rna = rseqt , meth = "equalfreq"){
  
  mar.T <- t(mar)
  rna.T <- t(rna)
  nbin <- sqrt(nrow(mar.T))
  mutual_lst <- c()
  for (i in 1:ncol(mar.T)) {
    mar.bin <- discretize(mar.T[ ,i], disc = meth, nbins = nbin)
    rna.bin <- discretize(rna.T[ ,i], disc = meth, nbins = nbin)
    info_i<- mutinformation(mar.bin, rna.bin)
    mutual_lst <- append(mutual_lst, info_i)
  }
  return(mutual_lst)
}

# Calculate fev, linear regression, scatter plots and named the picture
linear.scatter <- function(mar, rna, probes, func.lin, folderpath, num.pic, label){
  
  if (!file.exists(folderpath)) {
    dir.create(folderpath)
  } else {
    cat("The filepath already exists\n")
  }
  
  fev.lst <- c()
  stopifnot(num.pic < dim(mar)[1])
  if (missing(num.pic)) {num.pic = dim(mar)[1]} else { num.pic = num.pic}
  
  for (i in 1:(num.pic)) {
    # par(mfrow=c(1,2))
    fev <- lin.scatter(i, mar, rna, probes)
    fev.lst <- append(fev.lst, fev)
    if(label == TRUE){
      cat(sprintf("Press enter the label of %s functional probes (eg:Y/N): ", i))
      label_func <- readline()
      
      filename <- paste0(label_func, sprintf("_%s_Mutual_info = %s",i, round(probes$mutual[i], 3)),".png")
      file <- paste0('./',folderpath,"/",filename)
      png(file)
      scatter.lin(i, mar, rna, probes)
      dev.off()
    }}
  return(fev.lst)
}

# Scatter plot for typical probes
search.scatter <- function(mar, rna, search_gene, probes){
  for (i in search_gene){
    scatter(i, mar, rna, probes)
    readline("Please enter to continue...")
    }}
 
# Scatter plot
scatter <- function(i, mar, rna, probes){
  plot(mar[i, ], rna[i, ],
       xlab = "Microarray log expression",
       ylab = "RNA-seq log expression",
       xlim = marrt.lim,
       ylim = rseq.lim,
       main = sprintf("Gene = %s, Mutual_info = %.3f, R^2 = %.3f",
                      probes$genes[i], probes$mutual[i], probes$r.2[i]))
}

# Generate fev, prediction band based on linear regression and scatter plots
lin.scatter <- function(i, mar, rna, probes) {
  #par(mfrow=c(1,1))
  dev.new()
  dev.set(dev.list()[['RStudioGD']])
  scatter(i, mar, rna,  probes)
  
  lin.model <- lm(rna[i, ] ~ mar[i, ])
  abline(a = lin.model$coefficients[1], b = lin.model$coefficients[2], col = "red")

  fev.hat <- fev.func(rna[i, ], lin.model$fitted.values);
  
  idx <- order(mar[i, ]);
  pred.ym <- predict.lm(lin.model, interval = "confidence", level = 0.95)  # confidence interval for mean of y
  pred.yr <- predict.lm(lin.model, interval = "prediction", level = 0.95)  #
  lines(mar[i, ][idx], pred.ym[,"lwr"][idx], col = "blue", lty = "dashed")
  lines(mar[i, ][idx], pred.ym[,"upr"][idx], col = "blue", lty = "dashed")
  lines(mar[i, ][idx], pred.yr[,"lwr"][idx], col = "green", lty = "dashed")
  lines(mar[i, ][idx], pred.yr[,"upr"][idx], col = "green", lty = "dashed")
  
  readline("Press Enter to continue...")
  return(fev.hat)
}

#Generate fev, prediction band based on scam and scatter plots
curve.scatter <- function(mar, rna, probes, num.pic, meth = "mpi", k) {
  
  fev.lst <- c()
  if (missing(num.pic)) {num.pic = dim(mar)[1]} 
  else { 
    stopifnot(num.pic < dim(mar)[1])
    num.pic = num.pic}
  for (i in 1:(num.pic)) {
    dev.new()
    dev.set(dev.list()[['RStudioGD']])
    scatter(i, mar, rna,  probes)
    
    x1 <- mar[i, ];
    x2 <- rna[i, ];
    fit <- scam(x2 ~ s(x1, k = k, bs = meth));
    d <- data.frame(x1 = x1);
    x2.hat <- predict(fit, se.fit = TRUE, d);
    
    idx <- order(x1);
    lines(x1[idx], x2.hat$fit[idx], lwd = 2, type = "l",col = 'red')
    
    fev.hat <- fev.func(x2, fit$fitted.values);
    fev.lst <- append(fev.lst, fev.hat)
    readline("Press Enter to continue...")}
  return(fev.hat)}


# Plot scatter,fitted curves,confidence band for each curve probes
# Return FEV for each curve probes
#  fev_info <- function(mar , rna , gene, mutual, r, func = fev, meth = "mpi", alpha){
#   fev_lst <- c()
#   for (i in 1:1974) {
#     k <- 4;
#     x1 <- mar[i, ];
#     x2 <- rna[i, ];
#     fit <- scam(x2 ~ s(x1, k = k, bs = meth));
#     d <- data.frame(x1 = x1);
#     x2.hat <- predict(fit, se.fit = TRUE, d);
#     fev.hat <- func(x2, x2.hat$fit);
#     fev_lst <- append(fev_lst, fev.hat)
#     plot(x1, x2, xlim = marrt.lim, ylim = rseq.lim, 
#          main = sprintf("Gene = %s, Mutual_info = %.3f, R^2 = %.3f", 
#                         gene[i], mutual[i],r[i]))
#     # fitted curve
#     idx <- order(x1);
#     lines(x1[idx], x2.hat$fit[idx], lwd = 2, type = "l",col = 'red')
#     # confidence band
#     z <- x2 - x2.hat$fit
#     sd.pred <- (sum((z - mean(z))^2)/length(z))^(1/2)
#     a <- qnorm(1-alpha/2, mean = mean(z), sd = sd.pred)
#     upper <- x2.hat$fit[idx] + a
#     lower <- x2.hat$fit[idx] - a
#     lines(x1[idx], upper, pch=20, col = 'red' )
#     lines(x1[idx], lower, pch=20, col = 'red')
#     
#     print(paste("Processing item:", i))
#   #  Sys.sleep(1)
#   }
#   return(fev_lst)
# }
