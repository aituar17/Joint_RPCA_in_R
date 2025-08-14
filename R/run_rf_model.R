#CV evaluation

evaluate_model_cv <- function(features, labels, folds = 5, ntree = 500, seed = 42) {
    set.seed(seed)
    tab <- table(labels)
    if (length(labels) < 2L || length(tab) < 2L) stop("Need >=2 samples and >=2 classes.")
    folds <- max(2L, min(as.integer(folds), as.integer(min(tab)), length(labels) - 1L))
    folds_idx <- caret::createFolds(labels, k = folds, list = TRUE, returnTrain = FALSE)
  
    accs <- numeric(length(folds_idx)); aucs <- numeric(length(folds_idx))
    for (i in seq_along(folds_idx)) {
        test_idx  <- folds_idx[[i]]
        train_idx <- setdiff(seq_along(labels), test_idx)
        Xtr <- features[train_idx, , drop = FALSE]; Xte <- features[test_idx, , drop = FALSE]
        ytr <- labels[train_idx]; yte <- labels[test_idx]
        if (length(unique(ytr)) < 2L) { accs[i] <- NA_real_; aucs[i] <- NA_real_; next }
    
        pp <- prep_train_test(Xtr, Xte)
        rf <- randomForest(x = pp$Xtr, y = ytr, ntree = ntree)
    
        yhat <- predict(rf, pp$Xte, type = "response")
        accs[i] <- mean(yhat == yte)
    
        probs <- predict(rf, pp$Xte, type = "prob")
        all_lvls <- levels(labels)
        miss <- setdiff(all_lvls, colnames(probs))
        if (length(miss)) for (mm in miss) probs <- cbind(probs, setNames(rep(0, nrow(probs)), mm))
        probs <- probs[, all_lvls, drop = FALSE]
    
        if (length(unique(yte)) < 2L) {
            aucs[i] <- NA_real_
        } else {
            aucs[i] <- tryCatch(as.numeric(pROC::multiclass.roc(yte, probs)$auc), error = function(e) NA_real_)
        }
    }
    list(accuracy = mean(accs, na.rm = TRUE), auc = mean(aucs, na.rm = TRUE))
}