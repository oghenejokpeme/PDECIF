library(caret)
library(data.table)

estimate_regressor_performance <- function(y_actual, y_predicted, dec = 3){
  rmse      <- sqrt(sum((y_actual - y_predicted)^2) / length(y_actual))
  pearsonr  <- cor(y_actual, y_predicted, method = "pearson")
  error_estimates <- c(r = pearsonr, rmse = rmse)
}

get_ctrls <- function(){
  train_test <- trainControl(
    method = "LGOCV", 
    p = 0.7,
    number = 1,
    allowParallel = FALSE
  )

  cv5 <- trainControl(
    method = "cv",
    number = 5,
    allowParallel = FALSE
  )

  cv10 <- trainControl(
    method = "cv",
    number = 10,
    allowParallel = FALSE
  )

  list(traintest = train_test, cv5 = cv5, cv10 = cv10)
}

build_xgb_regressor <- function(x, y, ctrl){
  grid <- expand.grid(
    nrounds = c(500, 1000, 1500, 2000),
    max_depth = c(2, 4, 6, 8),
    eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
    gamma = 0.05,
    subsample = 0.5,
    colsample_bytree = 1,
    min_child_weight = 1
  )

  set.seed(47567)
  suppressWarnings(
    xgb_fit <- train(
      x = x, 
      y = y,
      method = "xgbTree",
      eval_metric = "rmse",
      verbose = 0, 
      trControl = ctrl,
      tuneGrid = grid,
      objective = "reg:squarederror"
    )
  )
  xgb_fit
}

read_response <- function(subset){
  basepath <- paste0("../input/y_", subset, "_")
  y_train <- read.csv(
    file = paste0(basepath, "train.csv"), 
    header = FALSE, 
    row.names = 1
  )
  y_train_vect <- y_train[, 2]
  names(y_train_vect) <- y_train[, 1]

  y_test  <- read.csv(
    file = paste0(basepath, "test.csv"), 
    header = FALSE, 
    row.names = 1
  )
  y_test_vect <- y_test[, 2]
  names(y_test_vect) <- y_test[, 1]

  list(train = y_train_vect, test = y_test_vect)
}

read_input <- function(subset, dataset, y){
  filepath <- paste0("../input/datasets/x_", dataset, "_", subset, ".csv")
  x <- read.csv(file = filepath, header = TRUE, row.names = 1)
  x_train <- x[names(y$train), ]
  x_test  <- x[names(y$test), ]

  list(train = x_train, test = x_test)
}

bootstrap_testing <- function(xgb_model, x_test, y_test, nrounds = 1000){
  set.seed(38854)
  total_samples <- nrow(x_test)
  results <- matrix(NA, nrow = nrounds, ncol = 2)
  for (i in 1:nrounds){
    bidx <- sample(1:total_samples, total_samples, replace = TRUE)
    bx_test <- x_test[bidx, ]
    by_test <- y_test[bidx]
    y_predicted <- predict(xgb_model, bx_test)
    results[i, ] <- estimate_regressor_performance(by_test, y_predicted)
  }
  colnames(results) <- c("R", "RMSE")
  results
}

# -----------------------------------------------------------------------------

main <- function(){
  datasets <- c(
    "ecif_with_hs",
    "spatial_ecif_with_hs"
  )
  subsets <- c("casf-07", "casf-13", "casf-16")
  subsets <- c("casf-16")
  ctrls <- get_ctrls()

  for (ctrl in names(ctrls)[1:1]){
    for (subset in subsets){
      results <- NULL
      y <- read_response(subset)
      for (dataset in datasets){
        note <- paste0(ctrl, " ", subset, " ", dataset)
        message(note)
        x <- read_input(subset, dataset, y)

        message("  building model")
        xgb_model <- build_xgb_regressor(x$train, y$train, ctrls[[ctrl]])
        saveRDS(
          xgb_model, 
          paste0("../output/models/", ctrl, "__", subset, "__", dataset, ".rds")
        )

        # Boostrap test performance.
        message("  bootstrap testing") 
        bootstrap_results <- bootstrap_testing(xgb_model, x$test, y$test, 50)
        write.csv(
          bootstrap_results, 
          file = paste0("../output/bootstrap_performance/", 
                        ctrl, "__", subset, "__", dataset, ".csv"),
          row.names = FALSE
        )  

        # General single test set performance.
        message("  general testing") 
        y_predicted <- predict(xgb_model, x$test)
        result <- estimate_regressor_performance(y$test, y_predicted)    
        results <- rbind(results, result)
      }
      rownames(results) <- datasets
      results <- round(results, 3)
      write.csv(
        results, 
        file = paste0("../output/general_performance/", ctrl, "__", subset, ".csv")
      )
    }
  }
}

main()