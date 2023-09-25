# Explore the matched expression data

library(io)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)


# microrray: Affymetrix U133 (fRMA normalized)
# RNA-seq: TCGA, featureCount

db <- org.Hs.eg.db;
source("./explore_func.R")
marr <- qread("micro_exprs_matched.rds");
rseq <- qread("rna_seq_exprs_t_matched.rds");


dim(marr) # 12664   294
dim(rseq) # 12664   294

# sanity check for matched rows and columns
stopifnot(rownames(marr) == rownames(rseq))
stopifnot(colnames(marr) == colnames(rseq))

# mapping gene names
entrezs <- rownames(marr);
genes <- mapIds(db, keys=entrezs, keytype="ENTREZID", column="SYMBOL");
names(genes) <- genes;

# check that there are no duplicated gene symbols
stopifnot(sum(duplicated(genes)) == 0)

rownames(marr) <- genes;
rownames(rseq) <- genes;


# Delete genes expressed 0 among all samples
uniq.num <- apply(rseq, 1, function(x) length(unique(x)))
uniq.lst <- which(uniq.num == 1)

marr <- marr[-uniq.lst, ]
rseq <- rseq[-uniq.lst, ]
genes <- genes[-uniq.lst]
entrezs <- entrezs[-uniq.lst];
dim(marr) # 12659  294
dim(rseq) # 12659   294


# Remove these genes which have high NA: pipeline deficiency?
rseq.p.na <- apply(rseq, 1, function(z) mean(is.na(z)));
hist(rseq.p.na, breaks=100)
na.idx <- rseq.p.na <= 0.5;
rseqt <- rseq[na.idx, ];
marrt <- marr[na.idx, ];
genes.f <- genes[na.idx];
entrezs.f <- entrezs[na.idx];
dim(marrt) # 12394   294
dim(rseqt) # 12394   294



marrt <- as.matrix(marrt);
rseqt <- as.matrix(rseqt);


smoothScatter(marrt, rseqt);
plot(marrt, rseqt, pch=".");


# Mean expressed for each probes
rseqt.m <- rowMeans(rseqt);
marrt.m <- rowMeans(marrt);
plot(marrt.m, rseqt.m, pch = ".");

# Range expressed for each probes
marrt.lim <- c(min(marrt), max(marrt));
rseq.lim <- c(min(rseqt, na.rm = TRUE), max(rseqt, na.rm = TRUE))



# ---

# Explore a way to measure the quality of probes 
# according to relation between Microarray and RNA-seq
# good probes have coherence between diff sequencing tools


# Correlation
marr.rseq.cors <- unlist(lapply(
  genes.f,
  function(gene) {
    cor(marrt[gene, ], rseqt[gene, ], use="complete.obs")
  }
));

hist(marr.rseq.cors, breaks=100) # Some Negative correlation 
mean(marr.rseq.cors^2 > 0.5, na.rm=TRUE)
head(sort(marr.rseq.cors, decreasing=TRUE))

# First explained variance by PCA
pca_lst <- pca.ratio(mar = marrt ,rna = rseqt)
hist(pca_lst, breaks = 100)


# Slope
slope_info <- unlist(lapply(genes.f,
  function(gene) {
    s_i <- cov(marrt[gene, ],rseqt[gene, ])/ var(marrt[gene, ])
  }))
hist(slope_info, breaks = 100 )


# Calculate Mutual information for each probes
mutual_I <- mutual_info(marrt, rseqt)
hist(mutual_I, breaks = 100 )


probes.info <- data.frame(
  genes = genes.f,
  entrezs = entrezs.f,
  marr_mean = marrt.m,
  rseq_mean = rseqt.m,
  corr = marr.rseq.cors,
  slope = slope_info,
  PCA_1 = pca_lst,
  mutual = mutual_I
  )

col_label <- c("low" = "blue", "high" = "red")
# Low correlation(less than 0.2) are mainly concentrated in the left 
ggplot(probes.info, aes(x = marr_mean, y = rseq_mean, colour = ifelse(corr < 0.2, "low", "high"))) +
  theme_classic() + geom_point(alpha=0.1) + scale_color_manual(values = col_label)

# Low pca_1(less than 0.8) are mainly concentrated in the upper left and lower left
ggplot(probes.info, aes(x = marr_mean, y = rseq_mean, colour =ifelse(PCA_1 < 0.8, "low", "high"))) +
  theme_classic() + geom_point(alpha=0.1) + scale_color_manual(values = col_label)

# Slope = 0 are more decentralized
ggplot(probes.info, aes(x = marr_mean, y = rseq_mean, colour = slope)) +
  theme_classic() + geom_point(alpha=0.1) + scale_colour_gradientn(colours = terrain.colors(10))

# Low mutual information(less than 0.2) are mainly concentrated in left
ggplot(probes.info, aes(x = marr_mean, y = rseq_mean, colour = ifelse(mutual_I < 0.1, "low", "high"))) +
  theme_classic() + geom_point(alpha=0.1) + scale_color_manual(values = col_label)

##Conclusion:
# 1. bad quality probes(non-functional probes) have low mean of Microarray intensify;


# Non-functional probes 
# CONDITION:  mutual_I < 0.05  
MI.cod <- quantile(mutual_I, 0.15)
nof.cond <- mutual_I < MI.cod
probes.nof <- subset(probes.info, nof.cond)
probes.nof <- subset(probes.info, nof.cond)
rseqt.nof <- rseqt[probes.nof$genes,  ]
marrt.nof <- marrt[probes.nof$genes,  ]
dim(probes.nof)[1] # 1859

# Functional probes
probes.f <- subset(probes.info, !nof.cond)
rseqt.f <- rseqt[rownames(rseqt) %in% probes.f$genes,  ]
marrt.f <- marrt[rownames(marrt) %in% probes.f$genes,  ]
dim(probes.f)[1] # 10535



# Linear regression
## Linear model
lm <- lapply(probes.f$genes, function(g) {lm(rseqt.f[g, ] ~ marrt.f[g, ])})
summary_lst <- lapply(lm, summary)

## linear property for each probes
rse <- r.2 <- adj.r.2 <- coef <- c()
for (i in 1:dim(probes.f)[1]){
  rse <- append(rse, summary_lst[[i]]$sigma) # Residual standard error
  r.2 <- append(r.2, summary_lst[[i]]$r.squared)
  adj.r.2 <- append(adj.r.2, summary_lst[[i]]$adj.r.squared)
}
probes.f <- cbind(probes.f, rse = rse, r.2 = r.2, adj.r.2 = adj.r.2)
hist(r.2, breaks = 100)


# Linear functional probes classified by R^2
# 237 probes follow linear model
curve.cond <- r.2 < 0.95
probes.f.lin <- subset(probes.f, !curve.cond)
curve_idx <- which(curve.cond)
rseqt.f.lin <- rseqt.f[-curve_idx, ]
marrt.f.lin <- marrt.f[-curve_idx, ]
func.lin <- lm[rownames(rseqt.f.lin)]
dim(rseqt.f.lin)[1]  # 237
length(func.lin) # 237
  


# ---
# Generate Scatter plots and fev for easy linear functional probes

# num.pic: draw num.pic pictures
# label: whether to named the picture(Y/N/S/HS/C)
# view: TRUE/FALSE (whether to view pic)
# press: TRUE/FALSE(whether see pic one by one )
# output fev.lin: Calculate fev 

fev.lin <- lin.scatter(marrt.f.lin, rseqt.f.lin, probes.f.lin, func.lin,
                 folderpath = "probes_plots", num.pic = 5, label = FALSE, 
                 level = 0.95, press = TRUE, view = TRUE)




# --- Explore nonlinear probes ---

# Easy curve functional probes
probes.f.cur <- subset(probes.f, curve.cond)
rseqt.f.cur <- rseqt.f[curve_idx, ]
marrt.f.cur <- marrt.f[curve_idx, ]
dim(marrt.f.cur) # 10298   294


# Classify curve probes into three categories:
## S curve / half-S curve / Inverted-C curve
plot(probes.f.cur$mutual, probes.f.cur$r.2)
r.2.cond <- quantile(probes.f.cur$r.2)[c(2,4)]
mutual.cond <- quantile(probes.f.cur$mutual)[c(2,4)]
abline(h = r.2.cond,  col = "red")
abline(v = mutual.cond, col = "green")


## Look at scatter of diff area according quantile of mutual and R^2
probes.h.m <- which(probes.f.cur$mutual > mutual.cond[2] & probes.f.cur$r.2 < r.2.cond[2]);
length(probes.h.m) # 417
search.scatter(marrt.f.cur, rseqt.f.cur, probes.h.m, probes.f.cur, pic_num = 10)


probes.l.h <- which(probes.f.cur$mutual < mutual.cond[1] & probes.f.cur$r.2 > r.2.cond[2])
length(probes.l.h) # 37
search.scatter(marrt.f.cur, rseqt.f.cur, probes.l.h, probes.f.cur, pic_num = 10)

probes.l.m <- which(probes.f.cur$mutual < mutual.cond[1] & probes.f.cur$r.2 < r.2.cond[2] & probes.f.cur$r.2 > r.2.cond[1] )
length(probes.l.m) # 370
search.scatter(marrt.f.cur, rseqt.f.cur, probes.l.m, probes.f.cur, pic_num = 10)

probes.m.l <- which(probes.f.cur$mutual < mutual.cond[2] & probes.f.cur$mutual > mutual.cond[1] & probes.f.cur$r.2 < r.2.cond[1])
length(probes.m.l) # 407
search.scatter(marrt.f.cur, rseqt.f.cur, probes.m.l, probes.f.cur, pic_num = 10)

probes.m.h <- which(probes.f.cur$mutual <  mutual.cond[2] & probes.f.cur$mutual > mutual.cond[1] & probes.f.cur$r.2 > r.2.cond[2])
length(probes.m.h) # 380
search.scatter(marrt.f.cur, rseqt.f.cur, probes.m.h, probes.f.cur, pic_num = 10)

probes.l.l <- which(probes.f.cur$mutual < mutual.cond[1] & probes.f.cur$r.2 < r.2.cond[1])
length(probes.l.l) # 2168
search.scatter(marrt.f.cur, rseqt.f.cur, probes.l.l, probes.f.cur, pic_num = 10)




# limited detection probes:
# "UBD", “ADAMDEC1", "NTS",  "EPYC", "ASS1",  "FABP6", "LYPD1"
#           .957      .952    .953     .954    .956     .971

# probe S-curve: signal saturation(rns exprs grow quick but micr grow slow)
# S-curve probes:
# "CYP4B1", "DEFB1"
# .973-S    .95-S



# -----

# Non-converge curve probes(.2%) in SCAM
non.conver.probes <- c(1975,2071,2127,2183,2796,2809,3421,4234,4882,5238,5247,5292,6124,6837,7601,7950,8176,8223,8359,8856,9176,9498) 
search.scatter(marrt.f.cur, rseqt.f.cur, non.conver.probes, probes.f.cur)
marrt.f.cur.noconv <- marrt.f.cur[non.conver.probes, ]
rseqt.f.cur.noconv <- rseqt.f.cur[non.conver.probes, ]
probes.f.cur.noconv <- probes.f.cur[non.conver.probes, ]
# rownames(marrt.f.cur)[non.conver.probes]
# noisy probes: "RASA2" "TRPV6"  "SCN8A" "KLRG1"  "CMKLR2" "APBB1IP" "BTRC" "RTEL1" "SLCO1A2" "FGF12"  "LTB4R" 
# non-noisy probes: "OPA3" "ITGAX" "RFWD3" "TRPV6"  "EFCC1"  "TMEM104" "RPL29" "ROCK2" "MED25" "WDFY3" "ZNF780B"

# Converge curve probes in SCAM
probes.f.cur.conv <- probes.f.cur[-non.conver.probes, ]
marrt.f.cur.conv <- marrt.f.cur[-non.conver.probes, ]
rseqt.f.cur.conv <- rseqt.f.cur[-non.conver.probes, ]

dim(marrt.f.cur.conv) # 10276   294


# Generate scatter and fev of curve probes based on SCAM 
# curve.scatter <- function(mar, rna, probes, num.pic, meth = "ts", k)
# num.pic: the nember of wanted printing picture;
# k: parameter of k in scam;

fev.curve <- curve.scatter(marrt.f.cur.conv, rseqt.f.cur.conv, 
                           num.pic = 10, probes.f.cur.conv, k = 4)


length(fev.curve)
fev <- fev.curve[1]
coef.scam <- fev.curve[2]


############TODO
# index of suspected heteroskedasticity probes by looking residual plots
# hetero <- c(8,263, 332, 342, 356, 360, 374,  397, 405, 406, 409, 435, 468, 482, 487, 510, 525, 551, 554, 555, 564)
# 8, 342,356, 360, 564
i = 564
Y = rseqt.f.cur.conv
X = marrt.f.cur.conv
y <- Y[i, ]
x <- X[i, ]
d <- data.frame(x = x)
fit1 <- scam(y ~ s(x, k = 20, bs = "mpi"));
pred1 <- predict(fit1, se.fit = TRUE,type="response", d)
scam.check(fit1)
k.check(fit1)


fit2 <- scam(y ~ s(x, k = 20, bs = "mpi"), weights =  residuals(fit1)^-2)
pred2 <- predict(fit2, se.fit = TRUE,type="response", d)
scam.check(fit2)
k.check(fit2)
AIC(fit1, fit2)


plot(x, y)
idx <- order(x)
lines(x[idx], pred1$fit[idx], col = "green")
lines(x[idx], pred2$fit[idx], col = "blue")


plot(x, residuals(fit1)/sd(residuals(fit1)), col = "green")
points(x, residuals(fit2)/sd(residuals(fit2)), col = "blue")


#############TODO
lis <-c (2,8,11,41,55)

residual.plot<- function(X, Y){
  for (i in 515:dim(X)[1]){
    y <- Y[i, ]
    x <- X[i, ]
    fit <- scam(y ~ s(x, k = 20, bs = "mpi"));
    # d <- data.frame(x = x)
    # pred <- predict(fit, d)
    idx <- order(x)
    #lines(x[idx], y.hat[idx])
    norm.res <- residuals(fit)/sd(residuals(fit))
    plot(x, norm.res, main = paste0("This is: ", i))
    abline(h = c(2, -2), col = 'red')
    print(paste0("This is: ", i))
    Sys.sleep(2)
  }}

residual.plot(Y = rseqt.f.cur.conv, X = marrt.f.cur.conv)



## conclusion
# 1. For linear probes, residual plot is not obvious, BP-test would be obvious;
# 2. There is a situation that weighted-scam reducing heteroskedasticity but AIC increased -- i = 564:
#      This is because of basis func and the number of it; It would changed if increased k;
# 3. Model fitting effect would be bad if weitht-scam but it does not have heteroskedasticity -- i = 564; 


## Solution:
# 1. filter sample with high and low quantile exp in Rna or Microarray would be better to reduce heteroskedasticity, because point is Centralized 
# 2. Quantitative methods to detect heteroskedasticity in nonlinear probes;
### TODO
# 1. Checking the effect of heteroskedasticity on scam models
##





# dev1.d <- fev.curve[[3]]
# dev1.se <- fev.curve[[4]]
# dev2.d <- fev.curve[[5]]
# dev2.se <- fev.curve[[6]]
# dev1.d.sum <- apply(dev1.d, 2, sum)
# dev1.d.mean <- apply(dev1.d, 2, mean)
# 
# library(Rtsne)
# 

# 
# coef.low.dim <- data.frame(
#   
#   x = coef.tsne$Y[ ,1],
#   y = coef.tsne$Y[ ,2]
# )
# 
# 
# ggplot(coef.low.dim, aes(x = x, y = y, color = dev1.d.mean)) +
#   geom_point()
# 
# hist(dev1.d.mean, breaks = 100)
# hist(dev1.d.sum, breaks = 100)
# length(which(dev1.d.mean > 2)) # 2.5%    249
# 
# sear.lst <- which(dev1.d.mean > 2)
# mar <- marrt.f.cur.conv[sear.lst, ]
# rna <- rseqt.f.cur.conv[sear.lst, ]
# probes <- probes.f.cur.conv[sear.lst, ]
#  
# for (i in 1:length(sear.lst)) {
#   fit <- scam(rna[i, ] ~ s(mar[i, ], k = 4, bs = "mpi"));
#   d <- data.frame(x1 = mar[i, ]);
#   plot(mar[i, ], rna[i, ])
#   Sys.sleep(1)
# }
# 



##########TODO
# non-converge probes



non.conver.probes
mar = marrt.f.cur.noconv
rna = rseqt.f.cur.noconv
probes = probes.f.cur.noconv

fev.lst <- c()
for (i in 1:dim(mar)[1]) {
  dev.new()
  dev.set(dev.list()[['RStudioGD']])
  
  i = 21
 
  x1 <- mar[i, ];
  x2 <- rna[i, ];
  fit <- scam(x2 ~ s(x1, bs = "mpi"), optimizer="bfgs");
  
  fit$aic
  fit$sp # estimated smoothing parameter
  fit$bfgs.info #
  
  
  plot(x1, fit$residuals)
  scam.check(fit)   # test model
  summary(fit)
  

  fit$scam
  fit$coefficients
  fit$aic
  fit$edf

  
  
  plot(fit) # check trend fo parameter 
  gam.check(fit)  # p-value of parameter to check the k is suitable
  
  d <- data.frame(x1 = x1);
  x2.hat <- predict(fit, se.fit = TRUE, d);
  
  scatter(i, mar, rna,  probes)
  idx <- order(x1);
  lines(x1[idx], x2.hat$fit[idx], lwd = 2, type = "l",col = 'red')
  
  fev.hat <- fev.func(x2, fit$fitted.values);
  fev.lst <- append(fev.lst, fev.hat)
  Sys.sleep(1)
  #  readline("Press Enter to continue...")}
}




### TODO
# classify probes category and fit on diff func
# Prediction interval on gam and scam


# ---

# TODO
# 0a. Find a way to reliably remove non-functional probes
# 0b. Classify probes into different classes based on response curves
#     (micoarray vs. RNA-seq log expression)
# 1. Find gene-specific maps for probes with linear response curves
# 2. Find gene-specific maps for probes with non-linear response curves
# 3a. Implement a method that takes a microarray log expression matrix and apply
#     the gene-specific maps in order to get the RNA-seq log expression matrix
# ----
# 3b. Apply this method to TCGA-OV microarray expression data, and you should get the
#     TCGA-OV RNA-seq expression data. Evaluate fit.
# 3c. Apply this method to other TCGA datasets with matched microarray and
#     RNA-seq and evaluate fit. If fit is not great, find out why and fix the
#     method.
# 3d. Apply this method to other expression datasets with matched microarray
#     and RNA-seq data (e.g. MAQC/SEQC)


































