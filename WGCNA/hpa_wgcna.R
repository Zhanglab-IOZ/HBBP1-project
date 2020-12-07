
library("WGCNA")

enableWGCNAThreads(nThreads=4)

hpa <- read.table("HPA_fpkm.txt", header=T, row.names=1)
sample_list <- read.table( "samplelist", header = T, row.names = 1)

sampleTree <- hclust(as.dist(1-cor(hpa)), method = "average")
pdf("SampleTree.pdf", height=10, width=10); plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex=0.5); dev.off()

hpa_filtered <- hpa[,!colnames(hpa) %in% c('ERR315467', 'ERR315375')] #ERR from sampling clusering outlier

shao=function(x){return(ifelse( sum(x)/length(x)>=0.5 , 1, 0)) }; t1 = apply( hpa_filtered, 1, function(x){sum( aggregate( ifelse(as.numeric(x)>=0.5,1,0) , by = list(a = sample_list[colnames(hpa_filtered),2] ), FUN = "shao" )[,2] ) >= 1 } )

hpa_filtered <- hpa_filtered[t1,]
#We required genes to be present in at least one tissue, with a mean FPKM across replicates higher than 0.5, which represented robust transcription.

tmp = apply( hpa_filtered , 1 , function(x){return(sd( as.numeric(x) ) / mean( as.numeric(x) ) ) } )

#quantile(tmp, seq(0,0.5,0.05))
#        0%         5%        10%        15%        20%        25%        30% 
#0.04747375 0.12069191 0.15234067 0.18283079 0.21986022 0.26995028 0.33555310 
#       35%        40%        45%        50% 
#0.42747807 0.54451842 0.68370812 0.84859238 

t2 = tmp >= 0.12
hpa_filtered <- hpa_filtered[t2,]
#We also set the coefficient of variation cutoff to 0.12 to remove genes less variable across tissues.

#dim(hpa_filtered)
#[1] 29122   169

write.table(hpa_filtered, "hpa_169sample-29122gene_fpkm.txt")

hpa_filtered <- log2(hpa_filtered +1)

powers <- c(1:20)
sft <- pickSoftThreshold(t(hpa_filtered), powerVector=powers, networkType="signed", verbose=3)

sft$powerEstimate
#[1] 7

netData <- blockwiseModules(datExpr=t(hpa_filtered), maxBlockSize=dim(hpa_filtered)[1], networkType="signed", power=sft$powerEstimate, mergeCutHeight=0.15, saveTOMFileBase = "hpa_169tissue_29122genes", saveTOMs=TRUE, minModuleSize=20, pamStage=FALSE, reassignThreshold = 1e-10, verbose=3, deepSplit=2)
#the core data of the network analysis

pdf("ModuleLabels.pdf", height=8, width=8); plotDendroAndColors(netData$dendrograms[[1]], netData$colors, "Module Colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05); dev.off()

geneTree= hclust(1-TOM,method="average")
tree <- cutreeHybrid(dendro = geneTree, pamStage=FALSE, minClusterSize = 150, cutHeight = 0.99999, deepSplit = 4, distM = as.matrix(1-TOM))
merged <- mergeCloseModules(exprData = t(hpa_filtered), colors = tree$labels, cutHeight = 0.15)
#cut the tree and merge some of the modules

mColorh <-labels2colors(merged$colors)

datTraits = matrix( 0 , nrow = dim( hpa_filtered )[2] , ncol = length( table(sample_list$SampleName.1 ) ) )
for( i in 1 : length( unique(sample_list$SampleName.1 ) ) ){ m = unique(sample_list$SampleName.1 )[i] ; datTraits[ which( sample_list[ colnames( hpa_filtered ) , 2 ] == m ) , i ] = 1 }
colnames( datTraits ) = unique(sample_list$SampleName.1 )
rownames( datTraits ) = colnames( hpa_filtered )
#rename the network matrix and add gene IDs and tissue IDs for rownames and colnames

#write.table(datTraits, "/rd/yuanh/hbbp1/wgcna/hpa/HPA_169samplematrix.txt", sep="\t")

newMEs_col <- sub('^..','',colnames(merged$newMEs))
merged_newMEs <- merged$newMEs[,paste('ME',as.character(sort(as.numeric(newMEs_col))),sep='')]
datTraits <- datTraits[, sort(colnames(datTraits))]

moduleTraitCor = cor( merged_newMEs , datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, dim(hpa_filtered)[2])
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");dim(textMatrix) = dim(moduleTraitCor)
#calculate the correlation between tissues and modules then the p value

png("hpa_Module-traitRelationships.png", width=1500, height=1200); labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(datTraits), yLabels = names( merged_newMEs ), ySymbols = names( merged_newMEs ), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1,1), main = paste("HPA Module-trait Relationships")); dev.off()

