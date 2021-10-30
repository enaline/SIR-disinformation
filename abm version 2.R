library(rsample)     # data splitting 
library(dplyr)       # data wrangling
library(rpart)       # performing regression trees
library(rpart.plot)  # plotting regression trees
library(ipred)       # bagging
library(caret)
library(MASS)
library(tree)
library(rattle)
library(RColorBrewer)
set.seed(1982)

basic <- read.csv("basic_out.csv", stringsAsFactors = F, header= T)
networked <- read.csv("network_degree.csv", stringsAsFactors = F, header= T)
soctype <- read.csv("soc_type_out.csv", stringsAsFactors = F, header= T)
spreaders <- read.csv("spreaders_last.csv", stringsAsFactors = F, header= T)


# split data (currently unnecessary) --------------------------------------
colnames(basic) <- c("X", "population_size", "maximum_threshold", "memory_length",
                     "infection_duration", "contact_frequency",
                     "initial_fraction", "red_ones", "max_step", "dose_unit")
basic <- basic[, -c(1, 9)]

basic_split <- initial_split(basic, prop = .7)
basic_train <- training(basic_split)
basic_test  <- testing(basic_split)



tb <- rpart(red_ones ~ ., data = basic_train, method = "anova")
rpart.plot(tb)
printcp(tb)
plotcp(tb)


tb2 <- predict(tb, data = basic_test)
RMSE(pred = tb2, obs = basic_train$red_ones)




###prune

pfit <- prune(tb, cp=tb$cptable[which.min(tb$cptable[,"xerror"]),"CP"])
pfit <- prune(tb, cp=0.01)
rpart.plot(pfit)


hyper_grid <- expand.grid(
  minsplit = seq(5, 20, 1),
  maxdepth = seq(8, 15, 1)
)

models <- list()

for (i in 1:nrow(hyper_grid)) {
  
  # get minsplit, maxdepth values at row i
  minsplit <- hyper_grid$minsplit[i]
  maxdepth <- hyper_grid$maxdepth[i]
  
  # train a model and store in the list
  models[[i]] <- rpart(
    red_ones ~ .,
    data    = basic,
    method  = "anova",
    control = list(minsplit = minsplit, maxdepth = maxdepth)
  )
}


get_cp <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  cp <- x$cptable[min, "CP"] 
}

get_min_error <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  xerror <- x$cptable[min, "xerror"] 
}


hyper_grid %>%
  mutate(
    cp    = purrr::map_dbl(models, get_cp),
    error = purrr::map_dbl(models, get_min_error)
  ) %>%
  arrange(error) %>%
  top_n(-5, wt = error)


opt.basic <- rpart(
  red_ones ~ .,
  data    = basic,
  method  = "anova",
  control = list(minsplit = 13,  cp = 0.01)
)
rpart.plot(opt.basic)




# no split basic----------------------------------------------------------------
head(basic)
dim(basic)
colnames(basic) <- c("X", "population_size", "maximum_threshold", "memory_length",
                     "infection_duration", "contact_frequency",
                     "initial_fraction", "red_ones", "max_step", "dose_unit")
basic <- basic[, -c(1, 9)]



m_basic<- rpart(red_ones ~ .,
                data    = basic,
                method  = "anova",
                control = c(list(cp = 0, xval = 10), rpart.control(minsplit = 40, minbucket = 20))
)

tiff("~/Desktop/csv's/g_basic.tiff", width = 4, height = 4, units = 'in', res = 720)
g_basic<-  rpart.plot(m_basic)
dev.off()


m_basic2<- rpart(red_ones ~ .,
                data    = basic,
                method  = "anova")


tiff("~/Desktop/csv's/g_basic.tiff", width = 4, height = 4, units = 'in', res = 720)
g_basic<-  rpart.plot(m_basic2)
dev.off()


plotcp(m_basic2)
abline(v = 12, lty = "dashed")
summary(m_basic2)

m_basic2$cptable[which.min(m_basic$cptable[,"xerror"]),"CP"]

hyper_grid <- expand.grid(
  minsplit = seq(5, 20, 1),
  maxdepth = seq(8, 15, 1)
)

models <- list()

for (i in 1:nrow(hyper_grid)) {
  
  # get minsplit, maxdepth values at row i
  minsplit <- hyper_grid$minsplit[i]
  maxdepth <- hyper_grid$maxdepth[i]
  
  # train a model and store in the list
  models[[i]] <- rpart(
    red_ones ~ .,
    data    = basic,
    method  = "anova",
    control = list(minsplit = minsplit, maxdepth = maxdepth)
  )
}


get_cp <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  cp <- x$cptable[min, "CP"] 
}

get_min_error <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  xerror <- x$cptable[min, "xerror"] 
}


hyper_grid %>%
  mutate(
    cp    = purrr::map_dbl(models, get_cp),
    error = purrr::map_dbl(models, get_min_error)
  ) %>%
  arrange(error) %>%
  top_n(-5, wt = error)


opt.basic <- rpart(
    red_ones ~ .,
    data    = basic,
  method  = "anova",
  control = list(minsplit = 13,  cp = 0.01)
)


tiff("~/Desktop/csv's/opt_basic.tiff", width = 4, height = 4, units = 'in', res = 720)
opt_basic<-  rpart.plot(opt.basic)
dev.off()

sort(opt.basic$variable.importance, decreasing =T)




# networked ---------------------------------------------------------------

head(networked)
dim(networked)


colnames(networked) <- c("X", "infection_duration","population_size","memory_length",
                         "contact_frequency", "maximum_threshold","proportion_of_spreaders",
                         "dose_unit", "society_type",  "average_degree", 
                          "red_ones")
networked <- networked[, -1]
networked$society_type = as.factor(networked$society_type) 


m_networked<- rpart(red_ones ~ .,
                    method="anova",
                    data= networked
)
rpart.plot(m_networked)

m_networked$cptable
printcp(m_networked)
plotcp(m_networked)
abline(v = 12, lty = "dashed")



##pruning
p_networked <- prune(m_networked, cp=m_networked$cptable[which.min(m_networked$cptable[,"xerror"]),"CP"])

rpart.plot(p_networked)
plotcp(p_networked)
abline(v = 12, lty = "dashed")



tiff("~/Desktop/csv's/g_pruned_networked.tiff", width = 4, height = 4, units = 'in', res = 720)
g_pruned_networked<-  rpart.plot(p_networked)
dev.off()


plotcp(p_networked)
summary(p_networked)

sort(p_networked$variable.importance, decreasing =T)





# hyper_networked ---------------------------------------------------------


hyper_grid <- expand.grid(
  minsplit = seq(5, 20, 1),
  maxdepth = seq(8, 15, 1)
)

models_n <- list()

for (i in 1:nrow(hyper_grid)) {
  
  # get minsplit, maxdepth values at row i
  minsplit <- hyper_grid$minsplit[i]
  maxdepth <- hyper_grid$maxdepth[i]
  
  # train a model and store in the list
  models_n[[i]] <- rpart(
    red_ones ~ .,
    data    = networked,
    method  = "anova",
    control = list(minsplit = minsplit, maxdepth = maxdepth)
  )
}


get_cp <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  cp <- x$cptable[min, "CP"] 
}

get_min_error <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  xerror <- x$cptable[min, "xerror"] 
}


hyper_grid %>%
  mutate(
    cp    = purrr::map_dbl(models_n, get_cp),
    error = purrr::map_dbl(models_n, get_min_error)
  ) %>%
  arrange(error) %>%
  top_n(-5, wt = error)


opt.networked <- rpart(
  red_ones ~ .,
  data    = networked,
  method  = "anova",
  control = list(maxdepth = 8, minsplit=6, cp = 0.01)
)


tiff("~/Desktop/csv's/opt_networked.tiff", width = 4, height = 4, units = 'in', res = 720)
opt_networked<-  rpart.plot(opt.networked)
dev.off()



# society type ------------------------------------------------------------

head(soctype)
dim(soctype)
colnames(soctype) <- c("X", "infection_duration", "population_size", "memory_length",
                       "contact_frequency","maximum_threshold", "proportion_of_spreaders",
                       "dose_unit", "society_type", "red_ones")
soctype <- soctype[, -1]
soctype$society_type <- as.factor(soctype$society_type)


m_soctype<- rpart(red_ones ~ .,
                  data    = soctype,
                  method="anova"
)

m_soctype$cptable

tiff("~/Desktop/csv's/g_soctype.tiff", width = 5, height = 5, units = 'in', res = 720)
g_soctype<-  rpart.plot(m_soctype)
dev.off()


plotcp(m_soctype)
summary(m_soctype)
sort(m_soctype$variable.importance, decreasing =T)


p_soc <- prune(m_soctype, cp=0.01, xerror =xerror )
rpart.plot(p_soc)

##pruning
p_soctype <- prune(m_soctype, cp=m_soctype$cptable[which.min(m_soctype$cptable[,"xerror"]),"CP"])
rpart.plot(p_soctype)
plotcp(p_soctype)
abline(v = 12, lty = "dashed")



tiff("~/Desktop/csv's/g_pruned_soctype.tiff", width = 4, height = 4, units = 'in', res = 720)
g_pruned_soctype<-  rpart.plot(p_soctype)
dev.off()





# spreaders ---------------------------------------------------------------


head(spreaders)
dim(spreaders)
colnames(spreaders) <- c("X", "infection_duration", "population_size", "memory_length",
                         "contact_frequency", "maximum_threshold", "proportion_of_spreaders",
                         "dose_unit", "red_ones")
spreaders <- spreaders[, -1]

m_spreaders<- rpart(red_ones ~ .,
                    data    = spreaders,
                    method  = "anova",
                    control = rpart.control(minsplit = 20, minbucket = 5)
)



tiff("~/Desktop/csv's/g_spreaders.tiff", width = 5, height = 5, units = 'in', res = 720)
g_spreaders<-  fancyRpartPlot(m_spreaders, sub ="", type= 1)
dev.off()

rpart.plot(m_spreaders)
plotcp(m_spreaders)
summary(m_spreaders)

sort(p_spreaders$variable.importance, decreasing =T)


##pruning
p_spreaders <- prune(m_spreaders, cp=m_spreaders$cptable[which.min(m_spreaders$cptable[,"xerror"]),"CP"])
rpart.plot(p_spreaders)
plotcp(p_spreaders)
summary(p_spreaders)
abline(v = 12, lty = "dashed")



tiff("~/Desktop/csv's/g_pruned_spreaders.tiff", width = 4, height = 4, units = 'in', res = 720)
g_pruned_spreaders<-  rpart.plot(p_spreaders)
dev.off()




# hyper_spreaders ---------------------------------------------------------


hyper_grid <- expand.grid(
  minsplit = seq(5, 20, 1),
  maxdepth = seq(8, 15, 1)
)

models_s <- list()

for (i in 1:nrow(hyper_grid)) {
  
  # get minsplit, maxdepth values at row i
  minsplit <- hyper_grid$minsplit[i]
  maxdepth <- hyper_grid$maxdepth[i]
  
  # train a model and store in the list
  models_s[[i]] <- rpart(
    red_ones ~ .,
    data    = spreaders,
    method  = "anova",
    control = list(minsplit = minsplit, maxdepth = maxdepth)
  )
}


get_cp <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  cp <- x$cptable[min, "CP"] 
}

get_min_error <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  xerror <- x$cptable[min, "xerror"] 
}


hyper_grid %>%
  mutate(
    cp    = purrr::map_dbl(models_s, get_cp),
    error = purrr::map_dbl(models_s, get_min_error)
  ) %>%
  arrange(error) %>%
  top_n(-5, wt = error)


opt.spreaders <- rpart(
  red_ones ~ .,
  data    = spreaders,
  method  = "anova",
  control = list(minsplit = 9,  cp = 0.01)
)


tiff("~/Desktop/csv's/opt_spreaders.tiff", width = 4, height = 4, units = 'in', res = 720)
opt_spreaders<-  rpart.plot(opt.spreaders)
dev.off()











# bagging -----------------------------------------------------------------

b_basic <- bagging(
  formula = red_ones ~ .,
  data    = basic,
  coob    = TRUE,
  nbagg=32
)

b_basic


# assess 10-50 bagged trees
ntree <- 10:50

# create empty vector to store OOB RMSE values
rmse <- vector(mode = "numeric", length = length(ntree))

for (i in seq_along(ntree)) {
  # reproducibility
  set.seed(1982)
  
  # perform bagged model
  model <- bagging(
    formula = red_ones ~ .,
    data    = basic,
    coob    = TRUE,
    nbagg   = ntree[i]
  )
  # get OOB error
  rmse[i] <- model$err
}

plot(ntree, rmse, type = 'l', lwd = 2)
abline(v = 25, col = "red", lty = "dashed")

