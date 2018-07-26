library(reshape2)
library(igraph)

source('networks_of_stranded_assets_library.R')

################################################################################
#
# Read in input-output data
#
################################################################################
io.dat <- read.csv("sample_fr_withhh.csv")
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
# #-------------------------------------------------------------------------------
# # Remove T & U
# #-------------------------------------------------------------------------------
# io.dat.mat <- remove_io(io.dat.mat, c("T","U"))
# 
# #-------------------------------------------------------------------------------
# # Remove P3 and P51G (have to add them first)
# #-------------------------------------------------------------------------------
# # First, add P3 and P51 as columns, equal to zero
# io.dat.mat[,c("P3","P51")] <- 0
# # Next, aggregate
# io.dat.mat <- agg_io(io.dat.mat, "FD", c("P3","P51"))

################################################################################
#
# Make A matrix and get rid of final demand
#
################################################################################
#-------------------------------------------------------------------------------
# Make the A & B matrices (including final demand)
#-------------------------------------------------------------------------------
io.z.mat <- t(io.dat.mat)
denom <- apply(io.z.mat, 1, sum)
io.a.mat <- t(t(io.z.mat)/denom) # Leontief
io.b.mat <- io.z.mat/denom # Ghosh

x.eig <- eigen(io.a.mat)
p.eig <- eigen(t(io.a.mat))
b.eig <- eigen(t(io.b.mat))

a.prim.vec <- data.frame(rownames(io.dat.mat),abs(as.numeric(x.eig$vectors[,1])),abs(as.numeric(p.eig$vectors[,1])),abs(as.numeric(b.eig$vectors[,1])))
names(a.prim.vec) <- c("sector", "x", "p","b")
a.prim.vec$net <- a.prim.vec$b - a.prim.vec$x
a.prim.vec$prod <- a.prim.vec$p * a.prim.vec$x
a.prim.vec[order(-a.prim.vec$x),]
a.prim.vec[order(-a.prim.vec$p),]
a.prim.vec[order(-a.prim.vec$b),]
a.prim.vec[order(-a.prim.vec$net),]
a.prim.vec[order(a.prim.vec$prod),]

with(a.prim.vec,plot(x,p,pch=""))
with(a.prim.vec,text(x,p,sector,cex=0.75))

with(a.prim.vec,plot(net,prod,pch=""))
with(a.prim.vec,text(net,prod,sector,cex=0.75))

with(a.prim.vec,plot(b,prod,pch=""))
with(a.prim.vec,text(b,prod,sector,cex=0.75))

#-------------------------------------------------------------------------------
# Get rid of final demand for analysis
#-------------------------------------------------------------------------------
# io.dat.mat <- remove_io(io.dat.mat, "FD")
# io.a.mat <- remove_io(as.data.frame(io.a.mat), "FD")
# io.b.mat <- remove_io(as.data.frame(io.b.mat), "FD")

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
# fbl.a <- fb_links(io.a.mat)
# fbl.b <- fb_links(io.b.mat)
# m.bias <- fbl.b$forward - fbl.a$backward
# print(names(m.bias[order(m.bias, decreasing=T)]))
# print(names(m.bias[order(fbl.b$forward, decreasing=T)]))


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

