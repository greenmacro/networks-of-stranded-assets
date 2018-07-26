library(reshape2)
library(igraph)

source('networks_of_stranded_assets_library.R')

################################################################################
#
# Read in input-output data
#
################################################################################
io.dat <- read.csv("data_fr.csv")
io.dat$total <- as.double(io.dat$total)

# Make into a matrix
# Note: matrix entries are the transpose of what is expected for the A matrix,
#    but correct for the Ghosh matrix:
#    io.dat.mat[i,j] = payment by sector i to sector j
io.dat.mat <- dcast(io.dat, purchasing ~ purchased)
rownames(io.dat.mat) <- io.dat.mat[,1]
io.dat.mat[,1] <- NULL

################################################################################
#
# Make Ghosh & Leontief matrices, with and without HH, extract final demand, and clean up
#
################################################################################
fd.bioph <- apply(io.dat.mat[c("HH","INV"), colnames(io.dat.mat) != "VA_tot"],2,sum)
fd.conv <- fd.bioph[names(fd.bioph) != "HH"]
io.dat.mat <- io.dat.mat[rownames(io.dat.mat) != "INV",]

# Double-counting by including both total VA and wages, so subtract:
x <- apply(io.dat.mat, 1, sum) - io.dat.mat$HH
# Now that it's been used to make the X vector, remove value added
io.dat.mat <- remove_io(io.dat.mat, "VA_tot")

g.bioph.mat <- t(t(io.dat.mat)/x)
a.bioph.mat <- t(io.dat.mat/x)

# For conventional matrix, remove households
io.dat.mat.nohh <- remove_io(io.dat.mat, "HH")
g.conv.mat <- t(t(io.dat.mat.nohh)/x[names(x) != "HH"])
a.conv.mat <- t(io.dat.mat.nohh/x[names(x) != "HH"])

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
fl.conv <- fb_links(g.conv.mat)$column
fl.biophys <- fb_links(g.bioph.mat)$column
print(names(fl.conv[order(fl.conv, decreasing=T)]))
print(names(fl.biophys[order(fl.biophys, decreasing=T)]))

bl.conv <- fb_links(a.conv.mat)$column
bl.biophys <- fb_links(a.bioph.mat)$column
print(names(bl.conv[order(bl.conv, decreasing=F)]))
print(names(bl.biophys[order(bl.biophys, decreasing=F)]))

net.conv <- fl.conv - bl.conv
print(names(net.conv[order(net.conv, decreasing=T)]))

net.biophys <- fl.biophys - bl.biophys
print(names(net.biophys[order(net.biophys, decreasing=T)]))

fl.cella.conv <- as.data.frame(fb_links_cella(a.conv.mat, fd.conv))
fl.cella.conv$total <- fl.cella.conv$forward + fl.cella.conv$backward
fl.cella.conv$net <- fl.cella.conv$forward - fl.cella.conv$backward
print(rownames(fl.cella.conv[order(fl.cella.conv$forward, decreasing=T),]))
print(rownames(fl.cella.conv[order(fl.cella.conv$backward, decreasing=T),]))
print(rownames(fl.cella.conv[order(fl.cella.conv$total, decreasing=T),]))
print(rownames(fl.cella.conv[order(fl.cella.conv$net, decreasing=T),]))

plot(fl.conv, bl.conv, pch="")
text(fl.conv, bl.conv, names(fl.conv), cex=0.75)

plot(fl.biophys, bl.biophys, pch="")
text(fl.biophys, bl.biophys, names(fl.biophys), cex=0.75)

#-------------------------------------------------------------------------------
# Inverted pyramid plot with minimal fully-connected network
#-------------------------------------------------------------------------------
# Remove households for conventional diagram, keep them in for "biophysical" version
m.mfc <- mfc(io.dat.mat) # mfc(io.dat.mat) # mfc(remove_io(io.dat.mat, "HH"))

m.mfc.graph <- graph_from_adjacency_matrix(as.matrix(m.mfc), mode = "directed", weighted = TRUE)
m.mfc.graph.noloops <- simplify(m.mfc.graph)
m.inout <- data.frame(cbind(degree(m.mfc.graph.noloops, mode = "in"),degree(m.mfc.graph.noloops, mode = "out")))
names(m.inout) <- c("in", "out")

sector <- "C19"

# Plot shells
shells <- ego_shells(m.mfc.graph.noloops, sector, "in")
shndx <- vector(mode = "integer", length = gorder(m.mfc.graph.noloops))
for (i in 1:length(shells)) {
  shndx[as.integer(V(m.mfc.graph)[shells[[i]]])] <- i
}
# assign colors
colbar <- rainbow(length(shells))
# create layout
V(m.mfc.graph.noloops)$size <- 0
V(m.mfc.graph.noloops)$shape <- "none"
V(m.mfc.graph.noloops)$label.cex <-0.75
E(m.mfc.graph.noloops)$arrow.size <- 0.25
ll <- ego_layout(m.mfc.graph.noloops, sector, "in", jitter=0.025, inverted = T)
plot(ego_directed(m.mfc.graph.noloops, sector, "in"), layout=ll)

