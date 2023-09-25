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


# Calculate mutual information for each probes between Microarray and RNA-seq
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
lin.scatter <- function(mar, rna, probes, func.lin, folderpath, num.pic , label, level, view, press){
  
  # folder of picture
  if (!file.exists(folderpath)) {
    dir.create(folderpath)
  } else {
    cat("The filepath already exists\n")
  }
  
  # number pictures
  if (missing(num.pic)) {num.pic = dim(mar)[1]; press = FALSE} else { num.pic = num.pic;  press = TRUE}
  stopifnot(num.pic <= dim(mar)[1])
  
  fev.lst <- c()
  for (i in 1:(num.pic)) {
    # Fit liner model and Prediction: PI
    lin.pred <- linear(i, x = mar, y = rna, alpha = level)
    
    # fev
    fev.hat <- fev.func(rna[i, ], lin.pred$model$fitted.values);
    fev.lst <- append(fev.lst, fev.hat)

    # Scatter plot with PI
    scatter(i, mar, rna,  probes, pred = lin.pred$pred)
    Sys.sleep(2)
    scatter(i, mar, rna,  probes, pred = lin.pred$wpred)
    Sys.sleep(2)
    # Picture by picture.
    if (press == TRUE & label ==TRUE) {
      
      
      cat(sprintf("Press enter the label of %s functional probes (eg:Y/N): ", i))
      label_func <- readline()
      
      # linear probes folder
      if(label_func == "Y"){
        
        lin.path <- paste0("./",folderpath, "linear_probes")
        dir.create(lin.path, recursive = TRUE, showWarnings = FALSE)
        filename <- paste0(label_func, sprintf("_%s_Mutual_info = %s",i, round(probes$mutual[i], 3)),".png")
        file <- paste0('./',folderpath,"/",filename)
        png(file)
      }
    } 
    else if(press == TRUE & label ==FALSE ){
      readline("Press Enter to continue...")
    } else {print(paste("Processing item:", i))}
  }
  return(fev.lst)
}



# Only Scatter plot for special gene index and quantities
search.scatter <- function(x, y, search_i, probes, pic_num = NULL){
  if (is.null(pic_num)){
    for (i in search_i){
      scatter(i, x, y, probes)
      readline("Please enter to continue...")
    }
  } 
  else { 
    stopifnot(pic_num > 0 & pic_num <= length(search_i))
    for (j in 1:pic_num) {
      scatter(j, x, y, probes)
      readline("Please enter to continue...")
    }
  }
}

# Scatter plot with or without PI and predicted model according <pred>
scatter <- function(i, x, y, probes, pred = NULL){
  
  gene.exp <- data.frame(
    "Microarray" = x[i, ], 
    "RNA-seq" = y[i, ]
    )
  data <- gene.exp
  
  pic <- ggplot(data, aes(x = Microarray, y = RNA.seq)) + geom_point() + 
    labs(x = "Microarray log expression", y = "RNA-seq log expression", 
         title = sprintf("Gene = %s, Mutual_info = %.3f, R^2 = %.3f",
                         probes$genes[i], probes$mutual[i], probes$r.2[i]))+
    xlim(marrt.lim) + ylim(rseq.lim) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (is.null(pred)){ print(pic) } 
  else {
    data <- cbind(gene.exp, pred)
    pic <- pic + geom_ribbon(data = data,mapping = aes(ymin = lwr, ymax = upr, x = Microarray), 
                      alpha = 0.4, inherit.aes = FALSE, fill = "#fdb338") +
           geom_line(data = data,mapping = aes(y = fit, x = Microarray), 
                inherit.aes = FALSE, size = 1, colour = "#025196")
    print(pic) 
    }
}  
  

# Generating linear model, CI, PI
linear <- function(i, x, y, alpha){
  
  model <- lm(y[i, ] ~ x[i, ])
  new <- data.frame(x[i, ])
  weight = residuals(model)^-2
  weight.model <- lm(y[i, ] ~ x[i, ], weights = weight)
  
  # LSE prediction band 
  y.hat <- predict.lm(model, newdata = new, interval = "prediction", alpha)
  y.pred <- predict.lm(model, newdata = new, interval = "prediction", alpha)
  
  # WLSE prediction band
  wy.pred <- predict.lm(weight.model, newdata = new, weights = weight, interval = "prediction", alpha)
  wy.hat <- predict.lm(weight.model, newdata = new, weights = weight, interval = "prediction", alpha)

  
  pred <- data.frame(
    fit = y.hat$fit[ ,1],
    lwr = y.pred$fit[ ,2],
    upr = y.pred$fit[ ,3]
  )
  weight.pred <- data.frame(
    fit = wy.hat$fit[ ,1],
    lwr = wy.pred$fit[ ,2],
    upr = wy.pred$fit[ ,3]
  )
  
  
  output <- list()
  output$model = model
  output$wmodel = weight.model
  output$pred = pred
  output$wpred = weight.pred

  return(output)
}



#Generate fev, prediction band based on scam and scatter plots
curve.scatter <- function(mar, rna, probes, num.pic, meth = "ts", k) {
  
  fev.lst <- c()
  if (missing(num.pic)) {num.pic = dim(mar)[1]} 
  else { 
    stopifnot(num.pic < dim(mar)[1])
    num.pic <- num.pic
  }
  coef.data <- data.frame()
  # deriv1.d <- c()
  # deriv1.se <- c()
  # deriv2.d <- c()
  # deriv2.se <- c()
  # s.lst <- c()
  # ho.lst <- c()
  # ht.lst <- c()
  for (i in 1:(num.pic)) {
    # dev.new()
    # dev.set(dev.list()[['RStudioGD']])

    
    x1 <- mar[i, ];
    x2 <- rna[i, ];
    fit <- scam(x2 ~ s(x1, k = k, bs = "mpi"));

    # plot(fit)
    # scam.check(fit)
    # summary(fit)
    
# ###TODO
#     
#     deriv1 <- derivative.scam(fit, smooth.number=1, deriv=1)
#     data <- data.frame(x = x1, y = x2)
#     idx <- sort(data$x, index = TRUE)
#     a <- data.frame(x = idx$x, y = deriv1$d[idx$ix])
#     plot1 <- ggplot(a) + geom_line() +
#                     aes(x = x, y = y) +  labs(title = "first derivate")
#    
#     #sum(deriv1$d)
#     deriv1.d <-cbind(deriv1.d, deriv1$d)
#     deriv1.se <- cbind(deriv1.se, deriv1$se.d)
#     
#     deriv2 <- derivative.scam(fit,smooth.number=1,deriv=2)
#     xx <- sort(data$x,index=TRUE)
#     b <- data.frame(x = xx$x, y = deriv1$d[xx$ix])
#     # plot2 <- ggplot(b) + geom_line() +
#     #   aes(x = x, y = y) +  labs(title = "second derivate")
#     
#     deriv2.d <-cbind(deriv2.d, deriv2$d)
#     deriv2.se <- cbind(deriv2.se, deriv2$se.d)
#     
#     #sum(deriv2$d)
#     # grid.arrange(plot1, plot2,  nrow = 1, ncol = 2)
#     # picname <- paste0(i," .png")
#     # ggsave(picname, width = 20, height = 10, dpi = 300)
# ####TODO    
    
    coef.data <- rbind(coef.data,fit$coefficients)
    
    d <- data.frame(x1 = x1);
    x2.hat <- predict_interval_scam(fit, d);
    
    # # scatter and prediction model
    scatter(i, mar, rna, probes, x2.hat)
    
    fev.hat <- fev.func(x2, fit$fitted.values);
    fev.lst <- append(fev.lst, fev.hat)
    print(paste("Processing item:", i))}
    #probes.class <- readline("Press the type of scatter trend(s/ht/ho):...")}
    # if (probes.class == s) {s.lst = append(s.lst, i) } 
    # else if (probes.class == ht) { ht.lst = append(ht.lst, i)} 
    # else {ho.lst = append(ho.lst, i)}
  return(list(fev.lst, coef.data))
  
}

# Code snippet from Dr.shih of explained-variance and thanks to him
# Original code: https://github.com/djhshih/explained-variance.git
# predict mean with prediction interval for new data
predict_interval_scam <- function(fit, newdata, level=0.95) {
  y.hat <- predict(fit, se.fit=TRUE, type="response", newdata);
  
  # sigma2.hat is based on old data (1..N)
  nu <- length(fit$residuals) - length(coef(fit));
  sigma2.hat <- 1 / nu * sum( fit$residuals^2 );
  
  # \sqrt( \hat{ var[w] } ) = 
  #   \sqrt( \hat{ var[y_{N+1}] } + \hat{ var[\hat{y_{N+1}}] }
  sqrt_hat_var_w <- sqrt(sigma2.hat + y.hat$se.fit^2);
  
  tq <- qt(1 - (1 - level)/2, df=nu);
  se <- tq*sqrt_hat_var_w;
  pred <- data.frame(
    fit = y.hat$fit,
    se.fit = y.hat$se.fit,
    lwr = y.hat$fit - se,
    upr = y.hat$fit + se
  )
  return(pred)
  }