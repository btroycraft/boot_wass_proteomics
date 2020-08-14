#Load in Data
library(readxl)
data <- read_excel("C:/Users/jamsh/OneDrive/Desktop/RTG_proteomics_Local/cDIA_MXLSAKBH-Exp1-2-3-4_Gill_with Volcanoplots_SelectedSubset (2).xlsx")



#transpose data
t_data = t(data)
rownames(t_data) <- c(1:nrow(t_data))
colnames(t_data) <- t_data[1,]
f_data = t_data[-c(1,2),]

#add environments label for each row
label = c(rep('M',24), rep('LS', 24), rep('AK', 24), rep('BH', 24))
labeled_data = cbind(f_data,label)
labeled_data = data.frame(labeled_data)



#Function to find pairwise sets of proteins that are significant based on Tukey HSD test and Fold-Change

find.included.proteins.pairwise <- function(labeled_data, p.val.cutoff, fc.cutoff){
  
  pairs.included = list(AK.BH = c(), AK.LS = c(), AK.M = c(), BH.LS = c(), BH.M = c(), LS.M = c())
  insig_diff = c()
  
  
  
  for (i in 1:(ncol(labeled_data)-1)){
    
    #run the ANOVA
    
    the.model = lm(as.numeric(labeled_data[1:nrow(labeled_data),i]) ~ labeled_data$label)
    anova.test = anova(the.model)
    
    #calculate fold change
    env.means = aggregate(as.numeric(labeled_data[1:nrow(labeled_data),i]), by = list(labeled_data$label), mean)
    lg_AK.BH.fc = log(env.means$x[1] / env.means$x[2],base = 2)
    lg_AK.LS.fc = log(env.means$x[1] / env.means$x[3], base = 2)
    lg_AK.M.fc = log(env.means$x[1] / env.means$x[4],base = 2)
    lg_BH.LS.fc = log(env.means$x[2] / env.means$x[3],base = 2)
    lg_BH.M.fc = log(env.means$x[2] / env.means$x[4], base = 2)
    lg_LS.M.fc = log(env.means$x[3] / env.means$x[4], base = 2)
    
    #if significant then run Tukey HSD and accumulate proteins with significant p-values
    
    if (anova.test[1,5] < 0.01){
      a <- aov(the.model)
      tuk.test = TukeyHSD(a)
      if(tuk.test$`labeled_data$label`[1,4] < p.val.cutoff && abs(lg_AK.BH.fc) > fc.cutoff){
        pairs.included$AK.BH = c(pairs.included$AK.BH, colnames(labeled_data)[i])
        
      }
      if(tuk.test$`labeled_data$label`[2,4] < p.val.cutoff && abs(lg_AK.LS.fc) > fc.cutoff){
        pairs.included$AK.LS = c(pairs.included$AK.LS, colnames(labeled_data)[i])
      }
      if(tuk.test$`labeled_data$label`[3,4] < p.val.cutoff && abs(lg_AK.M.fc) > fc.cutoff){
        pairs.included$AK.M = c(pairs.included$AK.M, colnames(labeled_data)[i])
      }
      if(tuk.test$`labeled_data$label`[4,4] < p.val.cutoff && abs(lg_BH.LS.fc) > fc.cutoff){
        pairs.included$BH.LS = c(pairs.included$BH.LS, colnames(labeled_data)[i])
      }
      if(tuk.test$`labeled_data$label`[5,4] < p.val.cutoff && abs(lg_BH.M.fc) > fc.cutoff){
        pairs.included$BH.M = c(pairs.included$BH.M, colnames(labeled_data)[i])
        
      }
      if(tuk.test$`labeled_data$label`[6,4] < p.val.cutoff && abs(lg_LS.M.fc) > fc.cutoff){
        pairs.included$LS.M = c(pairs.included$LS.M, colnames(labeled_data)[i])
      }
      
    } else{ 
      #add to the list of proteins with no significance
      insig_diff = append(colnames(labeled_data)[i], insig_diff)
  
    }
  }
  return(list(pairs.included = pairs.included,insig_diff = insig_diff))
}

#Make Data into form for EnhancedVolcano Package
construct.vp.data <- function(labeled_data){
  lg2FC.AK.BH = c()
  lg2FC.AK.LS = c()
  lg2FC.AK.M = c()
  lg2FC.BH.LS = c()
  lg2FC.BH.M = c()
  lg2FC.LS.M = c()

  pval.AK.BH = c()
  pval.AK.LS = c()
  pval.AK.M = c()
  pval.BH.LS = c()
  pval.BH.M = c()
  pval.LS.M = c()

  for (i in 1:(ncol(labeled_data)-1)){
    #Calculate log2FC
    env.means = aggregate(as.numeric(labeled_data[1:nrow(labeled_data),i]), by = list(labeled_data$label), mean)
    lg_AK.BH.fc = log(env.means$x[1] / env.means$x[2],base = 2)
    lg_AK.LS.fc = log(env.means$x[1] / env.means$x[3], base = 2)
    lg_AK.M.fc = log(env.means$x[1] / env.means$x[4],base = 2)
    lg_BH.LS.fc = log(env.means$x[2] / env.means$x[3],base = 2)
    lg_BH.M.fc = log(env.means$x[2] / env.means$x[4], base = 2)
    lg_LS.M.fc = log(env.means$x[3] / env.means$x[4], base = 2)

    #Make vectors of log2FC for each pair of environments
    lg2FC.AK.BH = c(lg2FC.AK.BH,lg_AK.BH.fc)
    lg2FC.AK.LS = c(lg2FC.AK.LS,lg_AK.LS.fc)
    lg2FC.AK.M = c(lg2FC.AK.M,lg_AK.M.fc)
    lg2FC.BH.LS = c(lg2FC.BH.LS,lg_BH.LS.fc)
    lg2FC.BH.M = c(lg2FC.BH.M,lg_BH.M.fc)
    lg2FC.LS.M = c(lg2FC.LS.M,lg_LS.M.fc)

    the.model = lm(as.numeric(labeled_data[1:nrow(labeled_data),i]) ~ labeled_data$label)
    a <- aov(the.model)
    tuk.test = TukeyHSD(a)
    pval.AK.BH = c(pval.AK.BH, tuk.test$`labeled_data$label`[1,4])
    pval.AK.LS = c(pval.AK.LS, tuk.test$`labeled_data$label`[2,4])
    pval.AK.M = c(pval.AK.M, tuk.test$`labeled_data$label`[3,4])
    pval.BH.LS = c(pval.BH.LS, tuk.test$`labeled_data$label`[4,4])
    pval.BH.M = c(pval.BH.M, tuk.test$`labeled_data$label`[5,4])
    pval.LS.M = c(pval.LS.M, tuk.test$`labeled_data$label`[6,4])

  }
  AK.BH.vpdata = data.frame("point.labels" = colnames(labeled_data[1:1505]), "log2FC" = lg2FC.AK.BH, "p-val" = pval.AK.BH)
  AK.LS.vpdata = data.frame("point.labels" = colnames(labeled_data[1:1505]), "log2FC" = lg2FC.AK.LS, "p-val" = pval.AK.LS)
  AK.M.vpdata = data.frame("point.labels" = colnames(labeled_data[1:1505]), "log2FC" = lg2FC.AK.M, "p-val" = pval.AK.M)
  BH.LS.vpdata = data.frame("point.labels" = colnames(labeled_data[1:1505]), "log2FC" = lg2FC.BH.LS, "p-val" = pval.BH.LS)
  BH.M.vpdata = data.frame("point.labels" = colnames(labeled_data[1:1505]), "log2FC" = lg2FC.BH.M, "p-val" = pval.BH.M)
  LS.M.vpdata = data.frame("point.labels" = colnames(labeled_data[1:1505]), "log2FC" = lg2FC.LS.M, "p-val" = pval.LS.M)

  return(list(AK.BH.vpdata, AK.LS.vpdata, AK.M.vpdata, BH.LS.vpdata, BH.M.vpdata, LS.M.vpdata))
}

#Create the Volcano Plots

library(EnhancedVolcano)
library(gridExtra)
library(grid)

vp.data = construct.vp.data(labeled_data)



p1 <- EnhancedVolcano(vp.data[[1]],
                lab = colnames(labeled_data[1:1505]),
                x = 'log2FC',
                y = 'p.val',
                xlim = c(-3, 4),
                title = 'Alaska versus Bodega Harbor',
                FCcutoff = 1,
                pointSize = 3.0,
                colAlpha = 0.5)

p2 <- EnhancedVolcano(vp.data[[2]],
                lab = colnames(labeled_data[1:1505]),
                x = 'log2FC',
                y = 'p.val',
                xlim = c(-3, 4),
                title = 'Alaska versus Lake Solano',
                FCcutoff = 1,
                pointSize = 3.0,
                colAlpha = 0.5)

p3 <- EnhancedVolcano(vp.data[[3]],
                lab = colnames(labeled_data[1:1505]),
                x = 'log2FC',
                y = 'p.val',
                xlim = c(-3, 4),
                title = 'Alaska versus Mexico',
                FCcutoff = 1,
                pointSize = 3.0,
                colAlpha = 0.5)

p4 <- EnhancedVolcano(vp.data[[4]],
                lab = colnames(labeled_data[1:1505]),
                x = 'log2FC',
                y = 'p.val',
                xlim = c(-3, 4),
                title = 'Bodega Harbor versus Lake Solano',
                FCcutoff = 1,
                pointSize = 3.0,
                colAlpha = 0.5)

p5 <- EnhancedVolcano(vp.data[[5]],
                lab = colnames(labeled_data[1:1505]),
                x = 'log2FC',
                y = 'p.val',
                xlim = c(-3, 4),
                title = 'Bodega Harbor versus Mexico',
                FCcutoff = 1,
                pointSize = 3.0,
                colAlpha = 0.5)

p6 <- EnhancedVolcano(vp.data[[6]],
                lab = colnames(labeled_data[1:1505]),
                x = 'log2FC',
                y = 'p.val',
                xlim = c(-3, 4),
                title = 'Lake Solano versus Mexico',
                FCcutoff = 1,
                pointSize = 3.0,
                colAlpha = 0.5)



grid.arrange(p1, p2,p3,p4,p5,p6,
             ncol=3,
             top = textGrob('EnhancedVolcano',
                            just = c('center'),
                            gp = gpar(fontsize = 1)))



#Intersections

#proteins which showed no significant differences among any of the environments
insig_prots = find.included.proteins.pairwise(labeled_data, p.val.cutoff = 0.01, fc.cutoff = 1)[2]


#Proteins that have a significant difference among some/all of the environments
sig_prots = find.included.proteins.pairwise(labeled_data, p.val.cutoff = 0.01, fc.cutoff = 1)

#Proteins that are significant particular to Alaska
AK.sig = Reduce(intersect, list(sig_prots[[1]]$AK.BH,sig_prots[[1]]$AK.LS, sig_prots[[1]]$AK.M))

#Proteins that are significant particular to Bodega Harbor
BH.sig = Reduce(intersect, list(sig_prots[[1]]$BH.LS,sig_prots[[1]]$AK.BH, sig_prots[[1]]$BH.M))

#Proteins that are significant particular to Lake Solano
LS.sig = Reduce(intersect, list(sig_prots[[1]]$BH.LS,sig_prots[[1]]$LS.M, sig_prots[[1]]$AK.LS))

#Proteins that are significant particular to Mexico
M.sig = Reduce(intersect, list(sig_prots[[1]]$LS.M,sig_prots[[1]]$AK.M, sig_prots[[1]]$BH.M))



#Take union of significant proteins
sig_prots_special = as.vector(Reduce(union, list(AK.sig, BH.sig, LS.sig, M.sig)))

#Plot Significant proteins in a different color
keyvals <- ifelse(
  vp.data[[1]]$point.labels %in% sig_prots_special, 'green',
         'black')
names(keyvals)[keyvals == 'black'] <- 'NotChosen'
names(keyvals)[keyvals == 'green'] <- 'Chosen'



p1_sig <- EnhancedVolcano(vp.data[[1]],
                      lab = colnames(labeled_data[1:1505]),
                      x = 'log2FC',
                      y = 'p.val',
                      xlim = c(-3, 4),
                      title = 'Alaska versus Bodega Harbor',
                      FCcutoff = 1,
                      pointSize = 3.0,
                      colAlpha = 0.5,
                      selectLab = rownames(vp.data[[1]])[which(names(keyvals) %in% c('Chosen'))],
                      colCustom = keyvals)





p2_sig <- EnhancedVolcano(vp.data[[2]],
                          lab = colnames(labeled_data[1:1505]),
                          x = 'log2FC',
                          y = 'p.val',
                          xlim = c(-3, 4),
                          title = 'Alaska versus Lake Solano',
                          FCcutoff = 1,
                          pointSize = 3.0,
                          colAlpha = 0.5,
                          selectLab = rownames(vp.data[[2]])[which(names(keyvals) %in% c('Chosen'))],
                          colCustom = keyvals)





p3_sig <- EnhancedVolcano(vp.data[[3]],
                          lab = colnames(labeled_data[1:1505]),
                          x = 'log2FC',
                          y = 'p.val',
                          xlim = c(-3, 4),
                          title = 'Alaska versus Mexico',
                          FCcutoff = 1,
                          pointSize = 3.0,
                          colAlpha = 0.5,
                          selectLab = rownames(vp.data[[3]])[which(names(keyvals) %in% c('Chosen'))],
                          colCustom = keyvals)





p4_sig <- EnhancedVolcano(vp.data[[4]],
                          lab = colnames(labeled_data[1:1505]),
                          x = 'log2FC',
                          y = 'p.val',
                          xlim = c(-3, 4),
                          title = 'Bodega Harbor versus Lake Solano',
                          FCcutoff = 1,
                          pointSize = 3.0,
                          colAlpha = 0.5,
                          selectLab = rownames(vp.data[[4]])[which(names(keyvals) %in% c('Chosen'))],
                          colCustom = keyvals)

p5_sig <- EnhancedVolcano(vp.data[[5]],
                          lab = colnames(labeled_data[1:1505]),
                          x = 'log2FC',
                          y = 'p.val',
                          xlim = c(-3, 4),
                          title = 'Bodega Harbor versus Mexico',
                          FCcutoff = 1,
                          pointSize = 3.0,
                          colAlpha = 0.5,
                          selectLab = rownames(vp.data[[5]])[which(names(keyvals) %in% c('Chosen'))],
                          colCustom = keyvals)

p6_sig <- EnhancedVolcano(vp.data[[6]],
                          lab = colnames(labeled_data[1:1505]),
                          x = 'log2FC',
                          y = 'p.val',
                          xlim = c(-3, 4),
                          title = 'Lake Solano versus Mexico',
                          FCcutoff = 1,
                          pointSize = 3.0,
                          colAlpha = 0.5,
                          selectLab = rownames(vp.data[[6]])[which(names(keyvals) %in% c('Chosen'))],
                          colCustom = keyvals)


grid.arrange(p1_sig, p2_sig,p3_sig,p4_sig,p5_sig,p6_sig,
             ncol=3,
             top = textGrob('EnhancedVolcano',
                            just = c('center'),
                            gp = gpar(fontsize = 1)))
