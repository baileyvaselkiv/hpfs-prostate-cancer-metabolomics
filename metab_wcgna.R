#working directory and library paths
getwd() 
.libPaths() 
options(stringsAsFactors=FALSE)

#load libraries
library(WGCNA)
library(survival)

#######################################################################################
#load metabolite data
#Variables: id, psatest, fastcat, matchid, chol_ab, fastingr, fastingm, fasting1, 
#fasting2, fasting3, timecatm, timecat2, timecat3, timecat4, yearcat2, yearcat3, 
#ageblood, time_bet, advcase, t3b_up, t4_up, dist_at_dx, bmi, lowgrade, higrade,
#higrd43, logrd34, 243 metabolites
metabData <- read.table("~/metabolomics/data/metab.for.ericka.2017.withgrade.csv", header=TRUE, sep=",")
colnames(metabData)

#load class information
classData <- read.table("~/metabolomics/data/metabolite_classes.2017.csv", header=TRUE, sep=",")
colnames(classData)

#load annotation information
annotData <- read.table("~/metabolomics/data/annotation.2017.csv", header=TRUE, sep=",")
colnames(annotData)

########################################################################################
#process metabolite data

#remove observations without matches
mids<-names(which(table(metabData$matchid)==2))
metabData<-metabData[metabData$matchid %in% mids,]
table(metabData$advcase)

#metabolite data
probeinds <- 5:247
data<-metabData[,probeinds]
names(data)

##################################################################################
#create expression and traits data for wcgna
metabExpr<-data
rownames(metabExpr)<-metabData$id
metabTraits<-metabData[,c("id","psatest", "fastcat", "matchid", "chol_ab", "fastingr", "fastingm", "fasting1", 
"fasting2", "fasting3", "timecatm", "timecat2", "timecat3", "timecat4", "yearcat2", "yearcat3", "ageblood", 
"time_bet", "advcase", "t3b_up", "t4_up", "dist_at_dx", "bmi", "lowgrade", "higrade", "higrd43", "logrd34")]

########################################################################################
#select thresholding power

#choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#call the network topology analysis function
sft = pickSoftThreshold(metabExpr, powerVector = powers, verbose = 5)
sft[[2]]

#plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
#scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
#this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#plot scale free topology fit
beta1=6
#the following computes the network connectivity (Connectivity) 
Connectivity= softConnectivity(metabExpr,power=beta1)
#create plot
par(mfrow=c(1,1)) 
scaleFreePlot(Connectivity, truncated=T,main= paste
              ("beta=",as.character(beta1))) 

########################################################################################
#construct network (detectCutHeight = 0.995)
net = blockwiseModules(metabExpr, maxBlockSize = 5000,
                       corType = "pearson",
                       power = 6, networkType = "unsigned", #unsigned adjacency = abs(cor)^power
                       TOMType = "signed", 
                       saveTOMs = FALSE, saveTOMFileBase = "~metabolomics/results/affyTOM",
                       deepSplit = 2, detectCutHeight = 0.99, minModuleSize = 20, #basic tree cut options
                       pamStage = FALSE, pamRespectsDendro = TRUE, #advanced tree cut options 
                       minKMEtostay = 0.3, minCoreKMESize = 20/3, minCoreKME = 0.5, #module trimming
                       reassignThreshold = 1e-6, #gene reassignment
                       mergeCutHeight = 0.15, #module merging
                       numericLabels = TRUE,
                       verbose = 3)
moduleLabels = net$colors
table(moduleLabels)
#convert labels to colors for plotting
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs
dendograms = net$dendrograms[[1]]
blockGenes = net$blockGenes[[1]]
#consider changing power, deepSplit, detectCutHeight, minModuleSize, mergeCutHeight

########################################################################################
#generate dendogram

#open a graphics window
sizeGrWindow(12, 9)
#plot the dendrogram and the module colors underneath
plotDendroAndColors(dendograms, moduleColors[blockGenes],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

########################################################################################
#save network
save(MEs, moduleLabels, moduleColors, dendograms, blockGenes, 
     file = "../erg_expression/results/affyNetworkConstruction-auto.RData")

########################################################################################
#visualize network

#load expression and phenotype data

#load network data
#lnames = load(file = "../erg_expression/results/affyNetworkConstruction-auto.RData")
#lnames

#select subset of genes for visualization
nMetabolites = ncol(metabExpr)
nSamples = nrow(metabExpr)
nSelect = 243
#for reproducibility, we set the random seed
set.seed(10)
select = sample(nMetabolites, size = nSelect)

#calculate dissimilary matrix and restrict to subset of genes
dissTOM=1-TOMsimilarityFromExpr(metabExpr, networkType = "unsigned", power=6, TOMType = "signed")
selectTOM = dissTOM[select, select]

#there's no simple way of restricting a clustering tree to a subset of genes, so we must re cluster
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

#open a graphical window
sizeGrWindow(9,9)
#taking the dissimilarity to a power makes the plot more informative by effectively changing 
#the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^10
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, subset of genes")

########################################################################################
#visualize relationship between clinical traits and modules

#recalculate module eigengenes
MEs = moduleEigengenes(metabExpr, moduleColors)$eigengenes

#calculate correlations and their p-values
moduleTraitCor = cor(MEs, metabTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#display correlations and their p-values within a heatmap plot
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metabTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

########################################################################################
#visualize relationship between modules

#recalculate module eigengenes
MEs = moduleEigengenes(metabExpr, moduleColors)$eigengenes

#isolate ERG status from the clinical traits
advcase = as.data.frame(metabTraits$advcase)
names(advcase) = "advcase"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, advcase))
MET = cbind(MEs, advcase)
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90, signed=TRUE, plotAdjacency=TRUE)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

#=====================================================================================
#
#  Code chunk 11: Identify genes with high significance and high module membership 
#
#=====================================================================================

#recalculate MEs with color labels
MEs0 = moduleEigengenes(metabExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, metabTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)




# Define variable weight containing the weight column of datTrait
weight = as.data.frame(metabTraits$advcase);
names(weight) = "erg_pos"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(metabExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(metabExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
module = "turquoise"
module = "blue"
column = match(module, modNames);
#moduleGenes = moduleColors==module;
moduleGenes = colorh3==module;


genelist<-as.matrix(colnames(metabExpr)[moduleGenes])

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for advanced case",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



pcaNames<-names(MEs)
x <- data.frame(MEs,
                id = metabData$id,
                psatest = metabData$psatest,
                fastcat = factor(metabData$fastcat, levels=c(1,2,3,4,NA), exclude=NULL),
                matchid = metabData$matchid,
                chol_ab = metabData$chol_ab,
                fastingr = metabData$fastingr,
                fastingm = metabData$fastingm,
                fasting1 = metabData$fasting1,
                fasting2 = metabData$fasting2,
                fasting3 = metabData$fasting3,
                timecatm = metabData$timecatm,
                timecat2 = metabData$timecat2,
                timecat3 = metabData$timecat3,
                timecat4 = metabData$timecat4,
                yearcat2 = metabData$yearcat2,
                yearcat3 = metabData$yearcat3,
                ageblood = metabData$ageblood,
                time_bet = metabData$time_bet,
                advcase = metabData$advcase,
                t3b_up = metabData$t3b_up,
                t4_up = metabData$t4_up,
                dist_at_dx = metabData$dist_at_dx,
                bmi = metabData$bmi,
                lowgrade = metabData$lowgrade,
                higrade = metabData$higrade,
                higrd43 = metabData$higrd43, 
                logrd34 = metabData$logrd34)

#main analysis
fits_clogit <- vector("list", length(pcaNames))
pvals_clogit <- rep(0, length(pcaNames))
estimates_clogit <- rep(0, length(pcaNames))

for (i in 1:length(pcaNames)) {
  #metabFormula<-as.formula(paste("advcase~",pcaNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",pcaNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=pcaNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:4,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_041417.csv",sep=',', col.names=NA)

#brown module info
brownresults<-as.matrix(colnames(metabExpr)[moduleColors=="brown"])
matchinds <- as.numeric(sapply(brownresults, function(x) which(annotData[,"name"] == x)))
brownresults <- cbind(brownresults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
