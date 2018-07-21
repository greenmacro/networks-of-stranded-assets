library(reshape2)
library(igraph)

source('networks_of_stranded_assets_library.R')

################################################################################
#
# Read in input-output data
#
################################################################################
io.dat <- read.csv("sample_fr.csv")
io.dat$total <- as.double(io.dat$total)

# Make into a matrix
# Note: matrix entries are the transpose of what is expected for the A matrix:
#    io.dat.mat[i,j] = payment by sector i to sector j
io.dat.mat <- dcast(io.dat, purchasing ~ purchased)
rownames(io.dat.mat) <- io.dat.mat[,1]
io.dat.mat[,1] <- NULL

################################################################################
#
# Clean up the data (this will depend on the data set & sector codes)
#
################################################################################
#-------------------------------------------------------------------------------
# Remove T & U
#-------------------------------------------------------------------------------
io.dat.mat <- remove_io(io.dat.mat, c("T","U"))

#-------------------------------------------------------------------------------
# Remove P3 and P51G (have to add them first)
#-------------------------------------------------------------------------------
# First, add P3 and P51 as columns, equal to zero
io.dat.mat[,c("P3","P51")] <- 0
# Next, aggregate
io.dat.mat <- agg_io(io.dat.mat, "FD", c("P3","P51"))

################################################################################
#
# Make A matrix and get rid of final demand
#
################################################################################
#-------------------------------------------------------------------------------
# Make the A matrix (including final demand)
#-------------------------------------------------------------------------------
io.a.mat <- t(io.dat.mat)
denom <- apply(io.a.mat, 1, sum)
io.a.mat <- io.a.mat/denom

#-------------------------------------------------------------------------------
# Get rid of final demand for analysis
#-------------------------------------------------------------------------------
io.dat.mat <- remove_io(io.dat.mat, "FD")
io.a.mat <- remove_io(as.data.frame(io.a.mat), "FD")

################################################################################
#
# Analysis
#
################################################################################
#-------------------------------------------------------------------------------
# Visualize matrix
#-------------------------------------------------------------------------------
plot_io(io.dat.mat)

#-------------------------------------------------------------------------------
# Backward and forward links
#-------------------------------------------------------------------------------
fbl <- fb_links(io.a.mat)
m.bias <- fbl$forward - fbl$backward
print(names(m.bias[order(m.bias, decreasing=T)]))


#-------------------------------------------------------------------------------
# Inverted pyramid plot with minimal fully-connected network
#-------------------------------------------------------------------------------
m.mfc <- mfc(io.dat.mat)

m.mfc.graph <- graph_from_adjacency_matrix(as.matrix(m.mfc), mode = "directed", weighted = TRUE)
m.mfc.graph.noloops <- simplify(m.mfc.graph)
m.inout <- data.frame(cbind(degree(m.mfc.graph.noloops, mode = "in"),degree(m.mfc.graph.noloops, mode = "out")))
names(m.inout) <- c("in", "out")

# Plot shells
shells <- ego_shells(m.mfc.graph.noloops, "B", "in")
shndx <- vector(mode = "integer", length = gorder(m.mfc.graph.noloops))
for (i in 1:length(shells)) {
  shndx[as.integer(V(m.mfc.graph)[shells[[i]]])] <- i
}
# assign colors
colbar <- rainbow(length(shells))
# create layout
V(m.mfc.graph.noloops)$size <- 0
V(m.mfc.graph.noloops)$label.cex <-1
E(m.mfc.graph.noloops)$arrow.size <- 0.25
ll <- ego_layout(m.mfc.graph.noloops, "B", "in", jitter=0.025, inverted = T)
plot(ego_directed(m.mfc.graph.noloops, "B", "in"), layout=ll)

