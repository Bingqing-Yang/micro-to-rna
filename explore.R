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


# Delete some genes expressed(RNA-seq) 0 among all samples
uni_count_idx <- c()
for (i in 1:nrow(rseq)){
  unique_count <- length(unique(rseq[i, ]))
  if (unique_count < 2){
    uni_count_idx <- append(uni_count_idx, i)
  }
}

marr <- marr[-uni_count_idx, ]
rseq <- rseq[-uni_count_idx, ]
genes <- genes[-uni_count_idx]
entrezs <- entrezs[-uni_count_idx];
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
# ---

# gene-specific (i.e. probe-specific) analysis

rseqt.m <- rowMeans(rseqt);
marrt.m <- rowMeans(marrt);
plot(marrt.m, rseqt.m, pch = ".");

marrt.lim <- c(min(marrt), max(marrt));
rseq.lim <- c(min(rseqt, na.rm = TRUE), max(rseqt, na.rm = TRUE))



# ---

# explore a way to automatically classify microarray probes

# Calculate correlation for each probes
marr.rseq.cors <- unlist(lapply(
  genes.f,
  function(gene) {
    cor(marrt[gene, ], rseqt[gene, ], use="complete.obs")
  }
));


hist(marr.rseq.cors, breaks=100)
mean(marr.rseq.cors^2 > 0.5, na.rm=TRUE)
head(sort(marr.rseq.cors, decreasing=TRUE))

# Calculate the ratio of first explained variance for each probes
pca_lst <- pca.ratio(mar = marrt ,rna = rseqt)
hist(pca_lst, breaks = 100)


# Calculate slope for each probes
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
);

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


# Non-functional probes 
# mutual_info of some probes show 0.     
mutual.nof.cond <- mutual_I < 0.05 
probes.nof <- subset(probes.info, mutual.nof.cond)
rseqt.nof <- rseqt[probes.nof$genes,  ]
marrt.nof <- marrt[probes.nof$genes,  ]
dim(probes.nof)[1] # 2160 

# Functional probes
probes.f <- subset(probes.info, !mutual.nof.cond)
rseqt.f <- rseqt[rownames(rseqt) %in% probes.f$genes,  ]
marrt.f <- marrt[rownames(marrt) %in% probes.f$genes,  ]
dim(probes.f)[1] # 10234



# Calculate linear regression for each probes
lm.fit <- lapply(probes.f$genes, function(g) {
 lm(rseqt.f[g, ] ~ marrt.f[g, ])
})
summary_lst <- lapply(lm.fit, summary)


rse <- r.2 <- adj.r.2 <- coef <- c()
for (i in 1:dim(probes.f)[1]){
  rse <- append(rse, summary_lst[[i]]$sigma) # Residual standard error
  r.2 <- append(r.2, summary_lst[[i]]$r.squared)
  adj.r.2 <- append(adj.r.2, summary_lst[[i]]$adj.r.squared)
  #  coef <- append(coef, summary_lst[[names(genes.fd[i])]]$coefficients)
}
probes.f <- cbind(probes.f, rse = rse, r.2 = r.2, adj.r.2 = adj.r.2)
hist(r.2, breaks = 100)



# Easy linear functional probes classified by R^2
r.2.curve.cond <- r.2 < 0.95
probes.f.lin <- subset(probes.f, !r.2.curve.cond)
curve_idx <- which(r.2.curve.cond)
rseqt.f.lin <- rseqt.f[-curve_idx, ]
marrt.f.lin <- marrt.f[-curve_idx, ]
func.lin <- lm.fit[rownames(rseqt.f.lin)]
dim(rseqt.f.lin)[1]  # 237
length(func.lin) # 237
  


# ---
# Generate Scatter plots and fev for easy linear functional probes

# num.pic: draw num.pic pictures
# label: named the picture(Y/N/S/HS/C)
# output fev.lin: Calculate fev 
fev.lin <- linear.scatter(marrt.f.lin, rseqt.f.lin, 
                 probes.f.lin, func.lin,
                 folderpath = "linear_probes_plots", 
                 num.pic = 10, label = FALSE)


# Easy curve functional probes
probes.f.cur <- subset(probes.f, r.2.curve.cond)
rseqt.f.cur <- rseqt.f[curve_idx, ]
marrt.f.cur <- marrt.f[curve_idx, ]
dim(marrt.f.cur) # 9997   294



# Deeper looking at curve probes
# Three categories: S curve / half-S curve / Inverted-C curve
plot(probes.f.cur$mutual, probes.f.cur$r.2)
r.2.cond <- quantile(probes.f.cur$r.2)[c(2,4)]
mutual.cond <- quantile(probes.f.cur$mutual)[c(2,4)]
abline(h = r.2.cond,  col = "red")
abline(v = mutual.cond, col = "green")


probes.h.m <- which(probes.f.cur$mutual > mutual.cond[2] & probes.f.cur$r.2 < r.2.cond[2]);
length(probes.h.m) # 411

# Same shape of scatter plot
probes.l.h <- which(probes.f.cur$mutual < mutual.cond[1] & probes.f.cur$r.2 > r.2.cond[2])
length(probes.l.h) # 35

probes.l.m <- which(probes.f.cur$mutual < mutual.cond[1] & probes.f.cur$r.2 < r.2.cond[2] & probes.f.cur$r.2 > r.2.cond[1] )
length(probes.l.m) # 360

probes.m.l <- which(probes.f.cur$mutual < mutual.cond[2] & probes.f.cur$mutual > mutual.cond[1] & probes.f.cur$r.2 < r.2.cond[1])
length(probes.m.l) # 395

probes.m.h <- which(probes.f.cur$mutual <  mutual.cond[2] & probes.f.cur$mutual > mutual.cond[1] & probes.f.cur$r.2 > r.2.cond[2])
length(probes.m.h) # 375

probes.l.l <- which(probes.f.cur$mutual < mutual.cond[1] & probes.f.cur$r.2 < r.2.cond[1])
length(probes.l.l) # 2104




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
# rownames(marrt.f.cur)[non.conver.probes]
# noisy probes: "RASA2" "TRPV6"  "SCN8A" "KLRG1"  "CMKLR2" "APBB1IP" "BTRC" "RTEL1" "SLCO1A2" "FGF12"  "LTB4R" 
# non-noisy probes: "OPA3" "ITGAX" "RFWD3" "TRPV6"  "EFCC1"  "TMEM104" "RPL29" "ROCK2" "MED25" "WDFY3" "ZNF780B"

# Converge curve probes in SCAM
probes.f.cur.conv <- probes.f.cur[-non.conver.probes, ]
marrt.f.cur.conv <- marrt.f.cur[-non.conver.probes, ]
rseqt.f.cur.conv <- rseqt.f.cur[-non.conver.probes, ]
dim(marrt.f.cur.conv)


# Generate scatter and fev of curve probes based on SCAM 
fev.curve <- curve.scatter(marrt.f.cur.conv, rseqt.f.cur.conv, 
                           probes.f.cur.conv, k = 4)
                          

# All curve probes is processed on GAM
mar = marrt.f.cur
rna = rseqt.f.cur
probes = probes.f.cur
fev.lst <- c()
for (i in 1:dim(marrt.f.cur)[1]) {
  dev.new()
  dev.set(dev.list()[['RStudioGD']])
  scatter(i, mar, rna,  probes)
  
  
  x1 <- mar[i, ];
  x2 <- rna[i, ];
  fit <- gam(x2 ~ s(x1));
  d <- data.frame(x1 = x1);
  x2.hat <- predict(fit, se.fit = TRUE, d);
  
  idx <- order(x1);
  lines(x1[idx], x2.hat$fit[idx], lwd = 2, type = "l",col = 'red')
  
  fev.hat <- fev.func(x2, fit$fitted.values);
  fev.lst <- append(fev.lst, fev.hat)
  Sys.sleep(1)
  #  readline("Press Enter to continue...")}
}



#### Summary:
# 1% curve probes do not converge on SCAM
# non-converge curve probe：some noisy and some non-noisy
# The predicted model on both scam and gam have good predictive effect
# GAM is not a robust model because k is hard to estimate, and diff k have diff model trend
# Some probes on GAM has decreasing trend

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


































