rm(list = ls())

# Load packages
library(sna)
library(igraph)
library(statnet)
library(PRROC) # To compute AUC (ROC) and  AUC (PR)

# Define working directory
setwd('~/Dropbox/PhD/2018 Summer/Network Analysis/hw1')

# Load data 
load('nigeria.rda')

# ------------ Part 1 ------------ #

# Check if sender can be the same as the receiver
nigeria[nigeria$sender==nigeria$receiver & nigeria$conflict > 0,]

# We have loops. To make the analysis simpler, we will remove the loops.

# Aggregate the data to sender and receiver
nig.small <- aggregate(nigeria$conflict, list(nigeria$sender, nigeria$receiver), sum, na.rm = TRUE)
names(nig.small) <-  c('sender', 'receiver', 'conflict')

# recode conflict = 1 if conflict > 0
nig.small$conflict <- ifelse(nig.small$conflict>0, 1, 0)

# Generate the matrix
nig.matrix <- as.matrix(table(nig.small)[, , 2])

# ------------ Part 2 ------------ #

# Given that we are dealing with conflict data, 
# I consider the most influential actor the one that had 
# the most direct participation in conflicts and is linked
# to actors that also had a large participation. Therefore,
# I use Eigenvector Centrality to measure the most
# the most influential actor.

# Import the matrix to igraph
g <- graph_from_adjacency_matrix(nig.matrix
	, mode = 'directed'
	, weighted = NULL
	, diag = FALSE)

# Calculate Eigenvector Centrality
ec <- eigen_centrality(g, directed = TRUE)$vector

# Get the maximum value
max(ec)

# Find all actors with EC equal to the maximum value
ec[which.max(ec)]

# Prepare data to run the block models
nig.matrix <- as_adjacency_matrix(g, type = 'both', edges = TRUE)
nig.matrix <- data.matrix(nig.matrix)
nigNet <- network(nig.matrix, directed = TRUE, loops = FALSE)

# Create the clusters
eclusters <- equiv.clust(nigNet)

# Generate a function to extract the block membership
block.membership <- function(x){
	m1 <- blockmodel(dat = nigNet, ec = eclusters, k = x)
	return(m1$block.membership[order(m1$order.vec)])
}


# Run block models to K = 2 to 13. Save the classification
bm.results <- sapply(2:14
	, block.membership
	, simplify = FALSE)

# We will use cross-validition to decide the parameter k.
# We will run logit models with the block membership as the unique 
# covariate

### Prepare the data
# Get the actors names
bmdf <- as.data.frame(cbind(rownames(nig.matrix)
	, Reduce('cbind', bm.results)))
names(bmdf) <- c('actor'
	, sapply(2:14, function(x) paste0('V', x)))
# # Merge the data
nig.small <- merge(x = nig.small
	, y = bmdf
	, by.x = 'sender'
	, by.y = 'actor')
nig.small <- merge(x = nig.small
	, y = bmdf
	, by.x = 'receiver'
	, by.y = 'actor')
# Generate sameblock dummies
genSameBlock <- function(n){
	ifelse(nig.small[, paste0('V', n, '.x')] == nig.small[, paste0('V', n, '.y')], 1, 0)	
}
sameBlock <- sapply(2:14, genSameBlock, simplify = FALSE)



# Function to run the cross validation
crossvalidation <- function(folds, foldn, blocklistn){

	# Start Cross-Validation
	tmp.df <- cbind(nig.small[, 'conflict'], sameBlock[[blocklistn]])
    train.df <- as.data.frame(tmp.df[-which(folds == foldn), ])
    test.df <- as.data.frame(tmp.df[which(folds == foldn), ])
    fit.model <- glm(train.df[, 1] ~ train.df[, 2], family = 'binomial')
    link <- coef(fit.model)[1] + (coef(fit.model)[2] * test.df[,2])
    pred <- 1/(1 + exp(-link))

	# Return a data.frame the predictions and original value
	return(cbind(test.df[,1], pred))
}

# Create folds
set.seed(1234)
folds <-  sample(rep(1:10, length.out = nrow(nig.small)))

# Run to each blockmodel results
cvRun <- sapply(1:13
	, function(x) sapply(1:10, crossvalidation, folds = folds, blocklistn = x)
	, simplify = FALSE)

# Combine the cross-validation runs
cvResult <- sapply(1:13, function(x) Reduce('rbind', cvRun[[x]]), simplify = FALSE)

# Generate a Function to get AUC (ROC) and  AUC (PR)
rocpr <- function(data){
	data <- as.data.frame(data)
	# Start to calculate the AUC (ROC) and  AUC (PR)
    fg <- data[data[, 1] == 1, 2]
	bg <- data[data[, 1] == 0, 2]
	# ROC     
	roc <- roc.curve(scores.class0 = fg, scores.class1 = bg)$auc
	# PR 
	pr <- pr.curve(scores.class0 = fg, scores.class1 = bg)$auc.integral
	# Return ROC and PR
	return(data.frame(roc, pr))
}

# Calculate the AUC (ROC) and  AUC (PR).
auc <- sapply(1:13, function(x) rocpr(cvResult[[x]])) 
auc

# Find the max
which.max(auc['roc',])
which.max(auc['pr',])

# We see that we should define k = 5. Now, we can produce the plot.

# Prepare vector with colors
block <- bmdf[, 'V5']

# Create group colors
colVec <- c('lightblue', 'pink', 'green', 'red', 'yellow')

# Assign colors to individual nodes based on block membership
bcols <- colVec[block]

# Now plot
set.seed(123)
plot(nigNet
	, displaylabels = TRUE
	, vertex.cex = 2
	, loop.cex = 4
    , label.cex = 0.8
    , edge.col = rgb(150, 150, 150, 100, maxColorValue = 255)
    , label.pos = 5
    , vertex.col = bcols)

# ------------ Part 3 ------------ #

# Add the block model group as a covariate:
network::set.vertex.attribute(nigNet
	, 'group'
	, as.numeric(as.character(block)))

# I hypothetize that group membership will
# decrease the likelihood of observing a edge.
# In other words, that actors in a same group will
# have a lower probability of fighting each other.
# I control for the number of egdes and popularity (idegree1.5)

# Now, we will run a ERGM.
m1 <- ergm(nigNet ~ edges
	+ idegree1.5
	+ nodematch('group'))
summary(m1)

# We can see that group is statistically significant and 
# has the expected sign. That is, when actors are in the same
# group, they have a lower probability of fighting each other.
# Actually, the probability of forming a edge creases in 0.17. 

# Now, we check the converge of the model using plots
mcmc.diagnostics(m1)

