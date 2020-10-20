setwd("C:/Users/jamsh/OneDrive/Desktop/RTG_proteomics_Local/boot_wass_proteomics-master/r")
library(TDA)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)


data <- read.table('cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv', header=TRUE, sep=',')

sub_inds <- read.csv("selected_subsets.csv", header = FALSE)

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


Names = rownames(data.frame(fit$points))

#Alaska MDS Plot
fit = isoMDS(D_ak, k = 2)
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



# fit = isoMDS(D_ak, k = 2)
# mds_ak_p <- plot(fit$points, main = "Alaska MDS Plot", type = "n")
# text(fit$points, labels = colnames(sub_prots), cex=.9)
# 
#fit = isoMDS(D_m, k = 2)
# mds_m_p <-plot(fit$points, main = "Mexico MDS Plot", type = "n")
# text(fit$points, labels = colnames(sub_prots), cex=.9)
# 
# fit = isoMDS(D_bh, k = 2)
# mds_bh_p <- plot(fit$points, main = "Bodega Harbor MDS Plot", type = "n")
# text(fit$points, labels = colnames(sub_prots), cex=.9)
#
#fit = isoMDS(D_ls, k = 2)
# mds_ls_p <-plot(fit$points, main = "Lake Solano MDS Plot", type = "n")
# text(fit$points, labels = colnames(sub_prots), cex=.9)

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
a <- data.frame(prots = ak_feature_prots, persistence.time = round(ripsd1.ak[, "Death"] - ripsd1.ak[, "Birth"], 4 ))
write.csv(a, "alaska_most_persistent_prots.csv")

b <- data.frame(prots = m_feature_prots, persistence.time = round(ripsd1.m[, "Death"] - ripsd1.m[, "Birth"], 4))
write.csv(b, "mexico_most_persistent_prots.csv")

c <- data.frame(prots = bh_feature_prots, persistence.time = round(ripsd1.bh[, "Death"] - ripsd1.bh[, "Birth"], 4))
write.csv(c, "bodega_most_persistent_prots.csv")

d <- data.frame(prots = ls_feature_prots, persistence.time = round(ripsd1.ls[, "Death"] - ripsd1.ls[, "Birth"] ))
write.csv(d, "lakesolano_most_persistent_prots.csv")

