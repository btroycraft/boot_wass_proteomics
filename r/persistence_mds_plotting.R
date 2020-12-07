setwd("C:/Users/jamsh/OneDrive/Desktop/RTG_proteomics_Local/boot_wass_proteomics-master/r")

library(TDA)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)
library(igraph)


data <- read.table('cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv', header=TRUE, sep=',')

sub_inds <- read.csv("../julia/selected_subsets.csv", header = FALSE)



ind_vec = as.vector(unlist(sub_inds[1])) + 2

sub_prots = cbind(data[ind_vec], location = data$location)


corrs_by_env = by(sub_prots[1:100], sub_prots$location, FUN = cor)


matrix_ones = matrix(rep(1, 100*100), 100, 100)


D_ak = matrix_ones - corrs_by_env$`Westchester Lagoon`
diag_ak = ripsDiag(D_ak, 1, 0.85, dist = "arbitrary")
a_p <- plot(diag_ak$diagram, main = "Alaska Population Persistence Diagram")
a_p2 <- plot(diag_ak$diagram, barcode = TRUE, main = "Alaska Population Persistence Diagram")




D_m = matrix_ones - corrs_by_env$`Laguna de la Bocana del Rosaria`
diag_m = ripsDiag(D_m, 1, 0.85, dist = "arbitrary")
m_p <- plot(diag_m$diagram, main = "Mexico Population Persistence Diagram")
m_p2 <- plot(diag_m$diagram, barcode = TRUE, main = "Mexico Population Persistence Diagram")




D_bh = matrix_ones - corrs_by_env$`Bodega Harbor`
diag_bh = ripsDiag(D_bh, 1, 0.85, dist = "arbitrary")
bh_p <- plot(diag_bh$diagram, main = "Bodega Harbor Population Persistence Diagram")
bh_p2 <- plot(diag_bh$diagram, barcode = TRUE, main = "Bodega Harbor Population Persistence Diagram")



D_ls = matrix_ones - corrs_by_env$`Lake Solano`
diag_ls = ripsDiag(D_ls, 1, 0.85, dist = "arbitrary")
ls_p <- plot(diag_ls$diagram, main = "Lake Solano Population Persistence Diagram")
ls_p2 <- plot(diag_ls$diagram, barcode = TRUE, main = "Lake Solano Population Persistence Diagram")




#Persistence Diagrams
par(mfrow=c(2,2))
plot(diag_ak$diagram, main = "Alaska Persistence Diagram")
plot(diag_m$diagram, main = "Mexico Persistence Diagram")
plot(diag_bh$diagram, main = "Bodega Harbor Persistence Diagram")
plot(diag_ls$diagram, main = "Lake Solano Persistence Diagram")


#Barcodes
# par(mfrow=c(2,2))
# plot(diag_ak$diagram, barcode = TRUE, main = "Alaska Barcodes")
# plot(diag_m$diagram, barcode = TRUE, main = "Mexico Barcodes")
# plot(diag_bh$diagram, barcode = TRUE, main = "Bodega Barcodes")
# plot(diag_ls$diagram, barcode = TRUE, main = "Lake Solano Barcodes")




#Alaska MDS Plot
fit = isoMDS(D_ak, k = 2)
Names = rownames(data.frame(fit$points))
mds_ak_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
  geom_point() +
  geom_text(aes(label=Names),hjust=0, vjust=.5) + geom_jitter() +
  ggtitle("Alaska MDS Plot") + xlim(-1,1.4)


#Mexico Plot
fit = isoMDS(D_m, k = 2)
mds_m_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
  geom_point() +
  geom_text(aes(label=Names),hjust=0, vjust=-.5) + geom_jitter() +
  ggtitle("Mexico MDS Plot") + xlim(-1,1.4)




#Bodega Plot
fit = isoMDS(D_bh, k = 2)
mds_bh_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
  geom_point() +
  geom_text(aes(label=Names),hjust=0, vjust=.5) + geom_jitter() +
  ggtitle("Bodega MDS Plot") + xlim(-1,1.2)


#Lake Solano Plot
fit = isoMDS(D_ls, k = 2)

mds_ls_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
  geom_point() +
  geom_text(aes(label=Names),hjust=0, vjust=.5) + geom_jitter() +
  ggtitle("Lake Solano MDS Plot") + xlim(-1,1.2)

g <- grid.arrange(mds_ak_p,mds_bh_p, mds_ls_p, mds_m_p,
                           ncol=2)





#Most persistent features

#Alaska most persistent
rips.ak <- ripsDiag(D_ak, 1,
                  0.85,
                  dist = "arbitrary",
                  library = "Dionysus",
                  location = TRUE)

ripsd1.ak <- rips.ak$diagram[(rips.ak$diagram)[, "dimension"] == 1, c("Birth", "Death")]
ord_ak <- order(ripsd1.ak[, "Death"] - ripsd1.ak[, "Birth"], decreasing = TRUE)
loc_ak = rips.ak$cycleLocation[(rips.ak$diagram)[, "dimension"] == 1]

ripsd1.ak <- ripsd1.ak[ord_ak, ]
loc_ak <- loc_ak[ord_ak]

#Mexico most persistent
rips.m <- ripsDiag(D_m, 1,
                   0.85,
                   dist = "arbitrary",
                   library = "Dionysus",
                   location = TRUE)

ripsd1.m <- rips.m$diagram[(rips.m$diagram)[, "dimension"] == 1, c("Birth", "Death")]
ord_m <- order(ripsd1.m[, "Death"] - ripsd1.m[, "Birth"], decreasing = TRUE)
loc_m = rips.m$cycleLocation[(rips.m$diagram)[, "dimension"] == 1]

ripsd1.m <- ripsd1.m[ord_m, ]
loc_m <- loc_m[ord_m]


#Bodega harbor most persistent
rips.bh <- ripsDiag(D_bh, 1,
                   0.85,
                   dist = "arbitrary",
                   library = "Dionysus",
                   location = TRUE)

ripsd1.bh <- rips.bh$diagram[(rips.bh$diagram)[, "dimension"] == 1, c("Birth", "Death")]
ord_bh <- order(ripsd1.bh[, "Death"] - ripsd1.bh[, "Birth"], decreasing = TRUE)
loc_bh = rips.bh$cycleLocation[(rips.bh$diagram)[, "dimension"] == 1]

ripsd1.bh <- ripsd1.bh[ord_bh, ]
loc_bh <- loc_m[ord_bh]


#Lake Solano most persistent
rips.ls <- ripsDiag(D_ls, 1,
                   0.85,
                   dist = "arbitrary",
                   library = "Dionysus",
                   location = TRUE)

ripsd1.ls <- rips.ls$diagram[(rips.ls$diagram)[, "dimension"] == 1, c("Birth", "Death")]
ord_ls <- order(ripsd1.ls[, "Death"] - ripsd1.ls[, "Birth"], decreasing = TRUE)
loc_ls = rips.ls$cycleLocation[(rips.ls$diagram)[, "dimension"] == 1]

ripsd1.ls <- ripsd1.ls[ord_ls, ]
loc_ls <- loc_ls[ord_ls]

#Get the protein names 
get_feature_prots <- function(loc){
  lapply(loc, function(x){
    prot_indices = unique(as.vector(x))
    sapply(prot_indices, function(inds){
      names(sub_prots[inds])
    })
    
  })
}

ak_feature_prots <- sapply(lapply(get_feature_prots(loc_ak), paste, collapse = "-"), as.vector)

m_feature_prots <- sapply(lapply(get_feature_prots(loc_m), paste, collapse = "-"), as.vector)

bh_feature_prots <- sapply(lapply(get_feature_prots(loc_bh), paste, collapse = "-"), as.vector)

ls_feature_prots <- sapply(lapply(get_feature_prots(loc_ls), paste, collapse = "-"), as.vector)


#Create Data frames and write to csv file
# a <- data.frame(prots = ak_feature_prots, persistence.time = round(ripsd1.ak[, "Death"] - ripsd1.ak[, "Birth"], 4 ))
# write.csv(a, "alaska_most_persistent_prots.csv")
# 
# b <- data.frame(prots = m_feature_prots, persistence.time = round(ripsd1.m[, "Death"] - ripsd1.m[, "Birth"], 4))
# write.csv(b, "mexico_most_persistent_prots.csv")
# 
# c <- data.frame(prots = bh_feature_prots, persistence.time = round(ripsd1.bh[, "Death"] - ripsd1.bh[, "Birth"], 4))
# write.csv(c, "bodega_most_persistent_prots.csv")
# 
# d <- data.frame(prots = ls_feature_prots, persistence.time = round(ripsd1.ls[, "Death"] - ripsd1.ls[, "Birth"], 4 ))
# write.csv(d, "lakesolano_most_persistent_prots.csv")
# 


#Plot the most persistent features by environment on an mds plot

#Get the most persistent groups of proteins
ak_mostpers <- get_feature_prots(loc_ak)[1:10]
bh_mostpers <- get_feature_prots(loc_bh)[1:10]
ls_mostpers <- get_feature_prots(loc_ls)[1:10]
m_mostpers <- get_feature_prots(loc_m)[1:10]

plot_mostpers_mds <- function(mostpers) {
  lapply(mostpers, function(prots){
    
    if (length(prots) == 0){
      NULL
    } else{
      a <- data.frame(cbind(sub_prots[prots], location = sub_prots$location))
      corrs = by(a[1:(length(a) - 1)], a$location, FUN = cor)
      n = dim(a)[2]
      matrix_ones = matrix(rep(1, (n-1)*(n-1)), n-1, n-1)
      dist_ak = matrix_ones - corrs$`Westchester Lagoon`
      dist_ls = matrix_ones - corrs$`Lake Solano`
      dist_m = matrix_ones - corrs$`Laguna de la Bocana del Rosaria`
      dist_bh = matrix_ones - corrs$`Bodega Harbor`
      
      
      fit = isoMDS(dist_ak, k = 2)
      Names = rownames(data.frame(fit$points))
      mds_ak_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
        geom_point() +
        geom_text(aes(label=Names), hjust = 0, vjust=.5) + geom_jitter() +
        ggtitle("Alaska MDS Plot") + xlim(-1,2.5) + ylim(-2.5, 2.5)
      
      fit = isoMDS(dist_bh, k = 2)
      mds_bh_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
        geom_point() +
        geom_text(aes(label=Names),hjust=0, vjust=.5) + geom_jitter() +
        ggtitle("Bodega Harbor MDS Plot") + xlim(-1,1.4)+ ylim(-2.5, 2.5)
      
      fit = isoMDS(dist_ls, k = 2)
      mds_ls_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
        geom_point() +
        geom_text(aes(label=Names),hjust=0, vjust=.5) + geom_jitter() +
        ggtitle("Lake Solano MDS Plot") + xlim(-1,1.4)+ ylim(-2.5, 2.5)
      
      fit = isoMDS(dist_m, k = 2)
      mds_m_p <- ggplot(data.frame(fit$points), aes(x = X1, y = X2, label = Names )) +
        geom_point() +
        geom_text(aes(label=Names),hjust=0, vjust=.5) + geom_jitter() +
        ggtitle("Mexico MDS Plot") + xlim(-2,3)+ ylim(-2.5, 2.5)
      g <- grid.arrange(mds_ak_p,mds_bh_p, mds_ls_p, mds_m_p,
                        ncol=2)
      return(corrs)

  }
})
}

vec_cor <- function(df){
  temp = df[1:ncol(df)]
  temp[lower.tri(temp)]
}

plot_mostpers_mds(ak_mostpers[1])
plot_mostpers_mds(bh_mostpers[7])
plot_mostpers_mds(ls_mostpers[1])
plot_mostpers_mds(m_mostpers[1])

#------------------------------PLOT NETWORKS FOR EACH OF THE PERSISTENT SUBSETS-------------------



##------------------------------------ALASKA MOST PERS------------------------------------------
#PLot MDS
a = plot_mostpers_mds(ak_mostpers[1])
comb.pairs = combn(ak_mostpers[[1]],2)
name.pairs = apply(comb.pairs, 2, paste, collapse="-")
edge_vec = unlist(strsplit(name.pairs, split = "-"))

par(mfrow = c(2,2))

#Create and plot network Alaska
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ak = vec_cor(data.frame(a[[1]]$`Westchester Lagoon`))
E(net)$width = abs(vec_corrs_ak) * 6
edge_colors = ifelse(abs(vec_corrs_ak) > 0.8, "red2", ifelse(abs(vec_corrs_ak) < 0.8 & abs(vec_corrs_ak) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Alaska Network")



#Create and plot network Bodega Harbor
net <- graph(edges = edge_vec, directed = F)
vec_corrs_bh = vec_cor(data.frame(a[[1]]$`Bodega Harbor`))
E(net)$width = abs(vec_corrs_bh) * 6
edge_colors = ifelse(abs(vec_corrs_bh) > 0.8, "red2", ifelse(abs(vec_corrs_bh) < 0.8 & abs(vec_corrs_bh) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Bodega Network")


#Create and plot network Lake Solano
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ls = vec_cor(data.frame(a[[1]]$`Lake Solano`))
E(net)$width = abs(vec_corrs_ls) * 6
edge_colors = ifelse(abs(vec_corrs_ls) > 0.8, "red2", ifelse(abs(vec_corrs_ls) < 0.8 & abs(vec_corrs_ls) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Lake Solano Network")


#Create and plot network Mexico
net <- graph(edges = edge_vec, directed = F)
vec_corrs_m = vec_cor(data.frame(a[[1]]$`Laguna de la Bocana del Rosaria`))
E(net)$width = abs(vec_corrs_m) * 6
edge_colors = ifelse(abs(vec_corrs_m) > 0.8, "red2", ifelse(abs(vec_corrs_m) < 0.8 & abs(vec_corrs_m) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Mexico Network")

dev.off()

##---------------------------------BODEGA MOST PERS----------------------------------------------
#PLot MDS
a = plot_mostpers_mds(bh_mostpers[7])
comb.pairs = combn(bh_mostpers[[7]],2)
name.pairs = apply(comb.pairs, 2, paste, collapse="-")
edge_vec = unlist(strsplit(name.pairs, split = "-"))

par(mfrow = c(2,2))

#Create and plot network Alaska
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ak = vec_cor(data.frame(a[[1]]$`Westchester Lagoon`))
E(net)$width = abs(vec_corrs_ak) * 6
edge_colors = ifelse(abs(vec_corrs_ak) > 0.8, "red2", ifelse(abs(vec_corrs_ak) < 0.8 & abs(vec_corrs_ak) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Alaska Network")



#Create and plot network Bodega Harbor
net <- graph(edges = edge_vec, directed = F)
vec_corrs_bh = vec_cor(data.frame(a[[1]]$`Bodega Harbor`))
E(net)$width = abs(vec_corrs_bh) * 6
edge_colors = ifelse(abs(vec_corrs_bh) > 0.8, "red2", ifelse(abs(vec_corrs_bh) < 0.8 & abs(vec_corrs_bh) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Bodega Network")


#Create and plot network Lake Solano
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ls = vec_cor(data.frame(a[[1]]$`Lake Solano`))
E(net)$width = abs(vec_corrs_ls) * 6
edge_colors = ifelse(abs(vec_corrs_ls) > 0.8, "red2", ifelse(abs(vec_corrs_ls) < 0.8 & abs(vec_corrs_ls) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Lake Solano Network")


#Create and plot network Mexico
net <- graph(edges = edge_vec, directed = F)
vec_corrs_m = vec_cor(data.frame(a[[1]]$`Laguna de la Bocana del Rosaria`))
E(net)$width = abs(vec_corrs_m) * 6
edge_colors = ifelse(abs(vec_corrs_m) > 0.8, "red2", ifelse(abs(vec_corrs_m) < 0.8 & abs(vec_corrs_m) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Mexico Network")


##--------------------------------Lake Solano MOST PERS----------------------------------------
#Find combinations of proteins
a = plot_mostpers_mds(ls_mostpers[1])
comb.pairs = combn(ls_mostpers[[1]],2)
name.pairs = apply(comb.pairs, 2, paste, collapse="-")
edge_vec = unlist(strsplit(name.pairs, split = "-"))

par(mfrow = c(2,2))

#Create and plot network Alaska
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ak = vec_cor(data.frame(a[[1]]$`Westchester Lagoon`))
E(net)$width = abs(vec_corrs_ak) * 6
edge_colors = ifelse(abs(vec_corrs_ak) > 0.8, "red2", ifelse(abs(vec_corrs_ak) < 0.8 & abs(vec_corrs_ak) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Alaska Network")



#Create and plot network Bodega Harbor
net <- graph(edges = edge_vec, directed = F)
vec_corrs_bh = vec_cor(data.frame(a[[1]]$`Bodega Harbor`))
E(net)$width = abs(vec_corrs_bh) * 6
edge_colors = ifelse(abs(vec_corrs_bh) > 0.8, "red2", ifelse(abs(vec_corrs_bh) < 0.8 & abs(vec_corrs_bh) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Bodega Network")


#Create and plot network Lake Solano
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ls = vec_cor(data.frame(a[[1]]$`Lake Solano`))
E(net)$width = abs(vec_corrs_ls) * 6
edge_colors = ifelse(abs(vec_corrs_ls) > 0.8, "red2", ifelse(abs(vec_corrs_ls) < 0.8 & abs(vec_corrs_ls) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Lake Solano Network")


#Create and plot network Mexico
net <- graph(edges = edge_vec, directed = F)
vec_corrs_m = vec_cor(data.frame(a[[1]]$`Laguna de la Bocana del Rosaria`))
E(net)$width = abs(vec_corrs_m) * 6
edge_colors = ifelse(abs(vec_corrs_m) > 0.8, "red2", ifelse(abs(vec_corrs_m) < 0.8 & abs(vec_corrs_m) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Mexico Network")


#----------------------------------------MEXICO MOST PERS---------------------------------------
#Find combinations of proteins
a = plot_mostpers_mds(m_mostpers[1])
comb.pairs = combn(m_mostpers[[1]],2)
name.pairs = apply(comb.pairs, 2, paste, collapse="-")
edge_vec = unlist(strsplit(name.pairs, split = "-"))

par(mfrow = c(2,2))

#Create and plot network Alaska
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ak = vec_cor(data.frame(a[[1]]$`Westchester Lagoon`))
E(net)$width = abs(vec_corrs_ak) * 6
edge_colors = ifelse(abs(vec_corrs_ak) > 0.8, "red2", ifelse(abs(vec_corrs_ak) < 0.8 & abs(vec_corrs_ak) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Alaska Network")



#Create and plot network Bodega Harbor
net <- graph(edges = edge_vec, directed = F)
vec_corrs_bh = vec_cor(data.frame(a[[1]]$`Bodega Harbor`))
E(net)$width = abs(vec_corrs_bh) * 6
edge_colors = ifelse(abs(vec_corrs_bh) > 0.8, "red2", ifelse(abs(vec_corrs_bh) < 0.8 & abs(vec_corrs_bh) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Bodega Network")


#Create and plot network Lake Solano
net <- graph(edges = edge_vec, directed = F)
vec_corrs_ls = vec_cor(data.frame(a[[1]]$`Lake Solano`))
E(net)$width = abs(vec_corrs_ls) * 6
edge_colors = ifelse(abs(vec_corrs_ls) > 0.8, "red2", ifelse(abs(vec_corrs_ls) < 0.8 & abs(vec_corrs_ls) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Lake Solano Network")


#Create and plot network Mexico
net <- graph(edges = edge_vec, directed = F)
vec_corrs_m = vec_cor(data.frame(a[[1]]$`Laguna de la Bocana del Rosaria`))
E(net)$width = abs(vec_corrs_m) * 6
edge_colors = ifelse(abs(vec_corrs_m) > 0.8, "red2", ifelse(abs(vec_corrs_m) < 0.8 & abs(vec_corrs_m) > 0.4, "orange2", "gray"))
coords = layout_in_circle(net, order = V(net))
plot(net, edge.color = edge_colors, layout = coords, main = "Mexico Network")



#-------------------------------------------------------------------------------------------------


##CLUSTER_PLOTTING

ak_clusters = read.csv("west (2).csv")
bh_clusters = read.csv("bodega.csv")
ls_clusters = read.csv("solano.csv")
m_clusters = read.csv("rosario.csv")


##Alaska Cluster MDS
#Get clustered proteins in a vector
ak_c_prots = ak_clusters$protein[which(ak_clusters$component == 1)]
bh_c_prots = bh_clusters$protein[which(bh_clusters$component == 1)]
ls_c_prots = ls_clusters$protein[which(ls_clusters$component == 1)]
m_c_prots = m_clusters$protein[which(m_clusters$component == 1)]


plot_mostpers_mds(list(ak_c_prots))
plot_mostpers_mds(list(bh_c_prots))
plot_mostpers_mds(list(ls_c_prots))
plot_mostpers_mds(list(m_c_prots))

