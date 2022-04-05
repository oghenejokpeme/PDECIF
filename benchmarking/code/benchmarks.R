library(caret)
library(data.table)

estimate_regressor_performance <- function(y_actual, y_predicted){
  rmse      <- sqrt(sum((y_actual - y_predicted)^2) / length(y_actual))
  pearsonr  <- cor(y_actual, y_predicted, method = "pearson")
  c(pearsonr, rmse)
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

read_response <- function(year){
  basepath <- paste0("../input/responses/y_", year, "_")
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

  if (length(intersect(names(y_train_vect), names(y_test_vect))) != 0){
    notification <- paste0("Shared samples between train and test in ",
                           year, "!")
    stop(notification)
  }

  list(train = y_train_vect, test = y_test_vect)
}

read_input <- function(year, dataset, y, dist = NULL){
  filepath <- ""
  if (is.null(dist)){
    filepath <- paste0("../input/datasets/", dataset, "_", year, ".csv")
  } else {
    filepath <- paste0("../input/datasets/", dataset, "_", year, "_", dist, ".csv")
  }
  
  x <- read.csv(file = filepath, header = TRUE, row.names = 1)
  x_train <- x[names(y$train), ]
  x_test  <- x[names(y$test), ]

  list(train = x_train, test = x_test)
}

# -----------------------------------------------------------------------------
# ECIF and PECIF independent experiments
general_experiments <- function(dtype, years, ctrls, dists = c(4, 6, 8, 10)){
  message(paste0("General experiments: ", dtype))

  for (year in years){
    message(paste0("  ", year))
    y <- read_response(year)

    for (dist in dists){
      x <- read_input(year, dtype, y, dist)

      for (ctrlname in names(ctrls)){
        message(paste0("    ", dist, " ", ctrlname))
        fname <- paste0(dtype, "_", year, "_", dist, "_", ctrlname)
        
        ctrl  <- ctrls[[ctrlname]]
        model <- build_xgb_regressor(x$train, y$train, ctrl)
        saveRDS(model, paste0("../output/models/", fname, ".rds"))

        y_predicted <- predict(model, x$test)
        performance <- estimate_regressor_performance(y$test, y_predicted)

        entry <- paste0(year, "-", dist, "\t", performance[1], "\t", performance[2])
        perfpath <- paste0("../output/results/", dtype, "_", ctrlname, ".txt")
        write(entry, append = T, file = perfpath) 
      }
    }
  }
}

# ECIF and PECIF + ligand features experiments
merge_datasets <- function(a, b, year){
  # Sanity check
  if (length(intersect(rownames(a$train), rownames(b$train))) != nrow(a$train)){
    notification <- paste0("Inconsistent samples in the merging of train sets ",
                           "in ", year)
    stop(notification)
  }
  if (length(intersect(rownames(a$test), rownames(b$test))) != nrow(a$test)){
    notification <- paste0("Inconsistent samples in the merging of test sets ",
                           "in ", year)
    stop(notification)
  }

  x_train <- cbind(a$train, b$train[rownames(a$train), ])
  x_test  <- cbind(a$test, b$test[rownames(a$test), ])

  list(train = x_train, test = x_test)
}

merged_experiments <- function(dtype, years, ctrls, dists = c(4, 6, 8, 10)){
  message(paste0("Merged experiments: ", dtype))

  for (year in years){
    message(paste0("  ", year))
    y <- read_response(year)
    ligand_x <- read_input(year, "ligand", y)

    for (dist in dists){
      dtype_x <- read_input(year, dtype, y, dist)
      x <- merge_datasets(ligand_x, dtype_x, year)

      for (ctrlname in names(ctrls)){
        message("    ", ctrlname)
        fname <- paste0("merged_ligand_", dtype, "_", year, "_", dist, "_", ctrlname)
        
        ctrl  <- ctrls[[ctrlname]]
        model <- build_xgb_regressor(x$train, y$train, ctrl)
        saveRDS(model, paste0("../output/models/", fname, ".rds"))

        y_predicted <- predict(model, x$test)
        performance <- estimate_regressor_performance(y$test, y_predicted)

        entry <- paste0(year, "-", dist, "\t", performance[1], "\t", performance[2])
        perfpath <- paste0("../output/results/merged_ligand_", dtype, "_", ctrlname, ".txt")
        write(entry, append = T, file = perfpath) 
      }
    }
  }
}

# -----------------------------------------------------------------------------
main <- function(){
  years <- c("casf-07", "casf-13", "casf-16", "casf-19")
  ctrls <- get_ctrls()

  general_experiments('ecif', years, ctrls)
  general_experiments('pecif', years, ctrls)
  merged_experiments('ecif', years, ctrls)
  merged_experiments('pecif', years, ctrls)
}

main()