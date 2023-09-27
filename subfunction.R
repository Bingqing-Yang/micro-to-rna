library(infotheo)
library(arules)
library(nlme)
library(mgcv)
library(scam)


#' Fev for Linear square and weighted linear square
#'
#' Calculate fev for linear square and weighted linear square;
#' Generate scatter plot with regression line for linear and weighted linear square
#'
#'
#' @param X  microarray data where rows are gene, columns are samples
#' @param Y  rna-seq data where rows are gene, columns are samples
#' @param probes  \code{data.frame}, probes information
#' @param folder folder for saving picture
#' @param num.pic number of output pictures
#' @param level confidence level
#' @param label whether to label pictures
#' @param view   whether to show pictures
#' @param press whether to press to show pictures
#' @param w weight of weighted linear model
#' @return list of fev values in linear model \code{fev.lin}, fev values in weighted linear model \code{fev.wlin}
#' @export
scatter.s <- function(X, Y, probes, pred, folder, num.pic, w, level, label, view, press) {
  # folder of picture
  if (!file.exists(folder)) {
    dir.create(folder)
  } else {
    cat("The filepath already exists\n")
  }

  # number pictures
  if (missing(num.pic)) {
    num.pic <- dim(X)[1]
    press <- FALSE
  } else {
    num.pic <- num.pic
    press <- TRUE
  }
  stopifnot(num.pic <= dim(X)[1])

  fev.lin <- c()
  for (i in 1:(num.pic)) {
    # # Fit liner model and Prediction: PI
    # lin.pred <- linear(i, X = X, Y = Y, level = level, weight = w)
    # 
    # # fev in model
    # fev.hat <- fev.func(Y[i, ], lin.pred$fit)
    # fev.lin <- append(fev.lin, fev.hat)
    # 
    # 
    # # Scatter plot with PI
    # predict.interval <- lin.pred[c(1:3)]
    predict.interval <- pred[[i]][c(1:3)]
    scatter(i, X, Y, probes, pred = predict.interval)


    # Picture by picture.
    if (press == TRUE & label == TRUE) {
      cat(sprintf("Press enter the label of %s functional probes (eg:Y/N): ", i))
      label_func <- readline()

      # folder for pictures
      if (label_func == "Y") {
        lin.path <- paste0("./", folder, "linear_probes")
        dir.create(lin.path, recursive = TRUE, showWarnings = FALSE)
        filename <- paste0(label_func, sprintf("_%s_Mutual_info = %s", i, round(probes$mutual[i], 3)), ".png")
        file <- paste0("./", folder, "/", filename)
        png(file)
      }
    } else if (press == TRUE & label == FALSE) {
      readline("Press Enter to continue...")
    } else {
      print(paste("Processing item:", i))
    }
  }
  # return(fev.lin)
}

# Generate fev, prediction band based on scam and scatter plots
curve.scatter <- function(mar, rna, probes, num.pic, meth = "ts", k) {
  fev.lst <- c()
  if (missing(num.pic)) {
    num.pic <- dim(mar)[1]
  } else {
    stopifnot(num.pic < dim(mar)[1])
    num.pic <- num.pic
  }
  coef.data <- data.frame()

  for (i in 1:(num.pic)) {

    x1 <- mar[i, ]
    x2 <- rna[i, ]
    fit <- scam(x2 ~ s(x1, k = k, bs = "mpi"))


    coef.data <- rbind(coef.data, fit$coefficients)

    d <- data.frame(x1 = x1)
    x2.hat <- predict_interval_scam(fit, d)

    # # scatter and prediction model
    scatter(i, mar, rna, probes, x2.hat)

    fev.hat <- fev.func(x2, fit$fitted.values)
    fev.lst <- append(fev.lst, fev.hat)
    print(paste("Processing item:", i))
  }
  return(list(fev.lst, coef.data))
}



#' Fraction of explained variance
#'
#' This function measure the fitting between true value \code{x}
#' and predicted value \code{x.hat}.
#'
#' 1 - sum(x - x.hat)^2/sum(x - mean(X))^2
#'
#' @param x         true value
#' @param x.hat     predicted value of model
#' @export
fev.func <- function(x, x.hat) {
  sse <- sum((x - x.hat)^2)
  sst <- sum((x - mean(x))^2)
  1 - sse / (sst + 1e-8)
}




#' Linear model or weighted linear model
#'
#' Implement linear model or weighted linear model according to the input parameter weight;
#' get prediction interval of specified model;
#'
#'
#' @param X  microarray data where rows are gene, columns are samples
#' @param Y  rna-seq data where rows are gene, columns are samples
#' @param i  i-th gene/probes
#' @param folder folder for saving picture
#' @param level confidence level
#' @param weight weight of model
#' @return list of models \code{model} and prediction interval \code{pred}
#' @export
linear <- function(i, X, Y, level, weight = NULL) {
  model <- lm(Y[i, ] ~ X[i, ], weights = weight)
  new <- data.frame(x = X[i, ])
  # weight <- residuals(model)^-2
  # weight.model <- lm(Y[i, ] ~ X[i, ], weights = weight)

  # prediction band
  if (is.null(weight)) {
    weight <- rep(1, ncol(X))
  }
  y.hat <- predict(model, newdata = new, weights = weight, interval = "prediction", level)
  y.pred <- predict(model, newdata = new, weights = weight, interval = "prediction", level)
  # # WLSE prediction band
  # wy.pred <- predict.lm(weight.model, newdata = new, weights = weight, interval = "prediction", level)
  # wy.hat <- predict.lm(weight.model, newdata = new, weights = weight, interval = "prediction", level)


  pred <- list(
    fit = y.hat$fit[, 1],
    lwr = y.pred$fit[, 2],
    upr = y.pred$fit[, 3],
    fit.se = y.hat$se.fit,
    model.res = model$residuals
  )
  # weight.pred <- data.frame(
  #   fit = wy.hat$fit[, 1],
  #   lwr = wy.pred$fit[, 2],
  #   upr = wy.pred$fit[, 3]
  # )


  # output <- list()
  # output$model <- model
  # # output$wmodel <- weight.model
  # output$pred <- pred
  # # output$wpred <- weight.pred


}

#' SCAM model
#'
#' Implement scam model or weighted scam model according whether input weight parameter
#' 
#'
#'
#' @param X  microarray data where rows are gene, columns are samples
#' @param Y  rna-seq data where rows are gene, columns are samples
#' @param i  i-th gene/probes
#' @param level confidence level
#' @param weight weight of model
#' @return list of response prediction \code{fit} and upper prediction interval \code{upr}
#'          and lower prediction interval \code{lwr}
#' @export
SCAM.model <- function(i, k, X, Y, level, weight = NULL) {

  samples.g <- data.frame(
    x = X[i, ],
    y = Y[i, ]
  )
  
  new <- data.frame(x = X[i, ])
  model <- scam(y ~ s(x, k = k, bs = "mpi"), data = samples.g) # tip: dim(model$model)[1] = dim(new)[1]
  y.pred <- predict_interval_scam(model, newdata = new, level)

  return(y.pred)
}


# # scatter and prediction model
# scatter(i, X, Y, probes, y.pred)


# Code snippet from Dr.shih of explained-variance and thanks to him
# Original code: https://github.com/djhshih/explained-variance.git
# predict mean with prediction interval for new data
predict_interval_scam <- function(model, newdata, level) {
  y.hat <- predict(model, se.fit = TRUE, type = "response", newdata)

  # sigma2.hat is based on old data (1..N)
  nu <- length(model$residuals) - length(coef(model))
  sigma2.hat <- 1 / nu * sum(model$residuals^2)

  # \sqrt( \hat{ var[w] } ) =
  #   \sqrt( \hat{ var[y_{N+1}] } + \hat{ var[\hat{y_{N+1}}] }
  sqrt_hat_var_w <- sqrt(sigma2.hat + y.hat$se.fit^2)

  tq <- qt(1 - (1 - level) / 2, df = nu)
  se <- tq * sqrt_hat_var_w
  pred <- data.frame(
    fit = y.hat$fit,
    se.fit = y.hat$se.fit,
    lwr = y.hat$fit - se,
    upr = y.hat$fit + se
  )
  return(pred)
}




#' Scatter plot for one gene
#'
#' Scatter plot given predicted result of model for one gene
#'
#'
#' @param X  microarray data where rows are gene, columns are samples
#' @param Y  rna-seq data where rows are gene, columns are samples
#' @param i  i-th gene/probes
#' @param probes  \code{data.frame}, probes information
#' @param pred \code{list}, predicted result of model containing fit, upper, lower
#' @param weight weight of model
#' @export
scatter <- function(i, X, Y, probes, pred) {
  gene.exp <- data.frame(
    "Microarray" = X[i, ],
    "RNA-seq" = Y[i, ]
  )
  data <- gene.exp

  pic <- ggplot(data, aes(x = Microarray, y = RNA.seq)) +
    geom_point() +
    labs(
      x = "Microarray log expression", y = "RNA-seq log expression",
      title = sprintf(
        "Gene = %s, Mutual_info = %.3f, R^2 = %.3f",
        probes$genes[i], probes$mutual[i], probes$r.2[i]
      )
    ) +
    xlim(marrt.lim) +
    ylim(rseq.lim) +
    theme(plot.title = element_text(hjust = 0.5))

  if (is.null(pred)) {
    print(pic)
  } else {
    data <- cbind(gene.exp, pred)
    pic <- pic + geom_ribbon(
      data = data, mapping = aes(ymin = lwr, ymax = upr, x = Microarray),
      alpha = 0.4, inherit.aes = FALSE, fill = "#fdb338"
    ) +
      geom_line(
        data = data, mapping = aes(y = fit, x = Microarray),
        inherit.aes = FALSE, size = 1, colour = "#025196"
      )
    print(pic)
  }
}

# Only Scatter plot for special gene index and quantities
search.scatter <- function(X, Y, search_i, probes, pred = NULL) {
  for (i in search_i) {
    scatter(i, X, Y, probes, pred)
    readline("Please enter to continue...")
  }
} 

#' Mutual information
#'
#' Calculate the mutual information between RNA and microarray for each probe
#'
#' @param X  Microarray data where rows are gene, columns are samples
#' @param Y  RNA-seq data where rows are gene, columns are samples
#' @return  \code{list} of mutual information for each genes
#' @export
mutual_info <- function(X, Y, meth = "equalfreq") {
  X.T <- t(X)
  Y.T <- t(Y)
  nbin <- sqrt(nrow(X.T))
  MI_lst <- c()
  for (i in 1:ncol(X.T)) {
    X.bin <- discretize(X.T[, i], disc = meth, nbins = nbin)
    Y.bin <- discretize(Y.T[, i], disc = meth, nbins = nbin)
    MI_i <- mutinformation(X.bin, Y.bin)
    MI_lst <- append(MI_lst, MI_i)
  }
  return(MI_lst)
}
