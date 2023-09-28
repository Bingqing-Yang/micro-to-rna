# Explore the matched expression data

library(io)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)


# microrray: Affymetrix U133 (fRMA normalized)
# RNA-seq: TCGA, featureCount

db <- org.Hs.eg.db
source("./subfunction.R")
marr <- qread("micro_exprs_matched.rds")
rseq <- qread("rna_seq_exprs_t_matched.rds")
dim(marr) # 12664   294
dim(rseq) # 12664   294

# sanity check for matched rows and columns
stopifnot(rownames(marr) == rownames(rseq))
stopifnot(colnames(marr) == colnames(rseq))

# mapping gene names
entrezs <- rownames(marr)
genes <- mapIds(db, keys = entrezs, keytype = "ENTREZID", column = "SYMBOL")
names(genes) <- genes
# check that there are no duplicated gene symbols
stopifnot(sum(duplicated(genes)) == 0)

rownames(marr) <- genes
rownames(rseq) <- genes
# Delete genes expressed 0 among all samples
uniq.num <- apply(rseq, 1, function(x) length(unique(x)))
uniq.lst <- which(uniq.num == 1)

marr <- marr[-uniq.lst, ]
rseq <- rseq[-uniq.lst, ]
genes <- genes[-uniq.lst]
entrezs <- entrezs[-uniq.lst]
dim(marr) # 12659  294
dim(rseq) # 12659   294


# Remove these genes which have high NA: pipeline deficiency?
rseq.p.na <- apply(rseq, 1, function(z) mean(is.na(z)))
hist(rseq.p.na, breaks = 100)
na.idx <- rseq.p.na <= 0.5
rseqt <- rseq[na.idx, ]
marrt <- marr[na.idx, ]
genes.f <- genes[na.idx]
entrezs.f <- entrezs[na.idx]
dim(marrt) # 12394   294
dim(rseqt) # 12394   294



marrt <- as.matrix(marrt)
rseqt <- as.matrix(rseqt)
smoothScatter(marrt, rseqt)
plot(marrt, rseqt, pch = ".")
# Mean expressed for each probes
rseqt.m <- rowMeans(rseqt)
marrt.m <- rowMeans(marrt)
plot(marrt.m, rseqt.m, pch = ".")
# Range expressed for each probes
marrt.lim <- c(min(marrt), max(marrt))
rseq.lim <- c(min(rseqt, na.rm = TRUE), max(rseqt, na.rm = TRUE))


# Calculate Mutual information for each probes
mutual_I <- mutual_info(marrt, rseqt)
hist(mutual_I, breaks = 100)



probes.info <- data.frame(
  genes = genes.f,
  entrezs = entrezs.f,
  marr_mean = marrt.m,
  rseq_mean = rseqt.m,
  mutual = mutual_I
)


# Non-functional probes
# CONDITION:  mutual_I < 0.05
MI.cod <- quantile(mutual_I, 0.15)
nof.cond <- mutual_I < MI.cod
probes.nof <- subset(probes.info, nof.cond)
probes.nof <- subset(probes.info, nof.cond)
rseqt.nof <- rseqt[probes.nof$genes, ]
marrt.nof <- marrt[probes.nof$genes, ]
dim(probes.nof)[1] # 1859



# Functional probes
probes.f <- subset(probes.info, !nof.cond)
rseqt.f <- rseqt[rownames(rseqt) %in% probes.f$genes, ]
marrt.f <- marrt[rownames(marrt) %in% probes.f$genes, ]
dim(probes.f)[1] # 10535



# Linear regression
## Linear model
lm <- lapply(probes.f$genes, function(g) {
  lm(rseqt.f[g, ] ~ marrt.f[g, ])
})
summary_lst <- lapply(lm, summary)

## linear property for each probes
rse <- r.2 <- adj.r.2 <- coef <- c()
for (i in 1:dim(probes.f)[1]) {
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
dim(rseqt.f.lin)[1] # 237
length(func.lin) # 237



# --- linear model ---

# Implement linear model for all genes
lin.pred <- lapply(1:dim(marrt.f.lin)[1], function(i) {
  linear(i, X = marrt.f.lin, Y = rseqt.f.lin, level = 0.95)
  }
)


# Calculate fev for all gene based on linear model
fev.lin <- unlist(lapply(1:length(lin.pred), function(i) {fev.func(rseqt.f.lin[i, ], lin.pred[[i]]$fit)}))


# Scatter plots of multiple genes fitted on the linear model
scatter.s(X = marrt.f.lin, Y = rseqt.f.lin, probes = probes.f.lin, pred = lin.pred,
               folder = "linear_plots", num.pic = 5, w = NULL,
               level = 0.95, label = FALSE, press = TRUE, view = TRUE
)
# num.pic: num.pic pictures drawn
# label: whether to named the picture(Y/N/S/HS/C)
# view: TRUE/FALSE (whether to view pic)
# press: TRUE/FALSE(whether see pic one by one )



# --- weighted least square ---

# Implement SCAM model for all genes
wlin.pred <- lapply(1:dim(marrt.f.lin)[1], function(i) {
  linear(i, X = marrt.f.lin, Y = rseqt.f.lin, level = 0.95, w = NULL)
}
)


# Calculate fev for all gene based on SCAM model
fev.lin <- unlist(lapply(1:length(wlin.pred), function(i) {fev.func(rseqt.f.lin[i, ], wlin.pred[[i]]$fit)}))


# Scatter plots of multiple genes fitted on the weighted linear model
scatter.s(X = marrt.f.lin, Y = rseqt.f.lin, probes = probes.f.lin, pred = wlin.pred,
               folder = "weighted_linear_plots", num.pic = 5, w = NULL,
               level = 0.95, label = FALSE, press = TRUE, view = TRUE
)



# --- SCAM model ---

# Nonlinear probes
probes.f.cur <- subset(probes.f, curve.cond)
rseqt.f.cur <- rseqt.f[curve_idx, ]
marrt.f.cur <- marrt.f[curve_idx, ]
dim(marrt.f.cur) # 10298   294


# Implement SCAM model for all genes
k.s <- c(seq(4, 20, 2))
scam.pred <- lapply(1:dim(marrt.f.cur)[1], function(i) {
  SCAM.model(i, k.s, X = marrt.f.cur, Y = rseqt.f.cur, level = 0.95)
  #print(i)
  }
)


# Calculate fev of scam model
fev.scam <- unlist(lapply(1:dim(marrt.f.cur)[1], function(i) 
  {fev.func(rseqt.f.cur[i, ], scam.pred[[i]]$fit)}
  )
)
hist(fev.scam, breaks = 100)



# Scatter plots of multiple genes fitted on SCAM model
scatter.s(X = marrt.f.cur, Y = rseqt.f.cur, probes = probes.f.cur, pred = scam.pred,
          folder = "scam_plots", num.pic = 5, w = NULL,
          level = 0.95, label = FALSE, press = TRUE, view = TRUE
)













# TODO
# 1. check the edf of non-conver probes, see does it equal to 1;


