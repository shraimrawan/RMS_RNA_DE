rm (list = ls())
library(limma)
library(edgeR)

data= read.table("genecount_09.07.20.txt", header= TRUE, sep ="\t")
head(data)
dim(data)
group= as.factor(c("ERMS","ARMS","ARMS","ERMS","ARMS","ARMS","ERMS","ARMS","ERMS","ERMS","ARMS"))

gene_anot=read.table('mart_export.txt', header=TRUE, sep=',')

normal_data=read.table('samplegenes2.gct', sep="\t", header=TRUE, check.names = FALSE)
head(normal_data)
tissue=read.table("phenotype_GTEX.txt", sep="\t", header=TRUE)
head(tissue)
samples=colnames(normal_data)
tissue_id=tissue$SAMPID
i=match(samples[], tissue_id)
i=i[-(1:2)]
tissue_order=tissue$SMTS[i]

##Filter Expression

g=match(normal_data$Name, data$geneid)
sample_data=data[g,]
all_data=merge(sample_data, normal_data, by.x = "geneid", by.y="Name", all.x = TRUE)
all_data=na.omit(all_data)

rownames(all_data)=all_data$geneid
all_data$Description= NULL
all_data$geneid=NULL
cpm=cpm(all_data)
lcpm = cpm(all_data, log=TRUE)
group_all=as.factor(c(as.character(group), as.character(tissue_order)))
group_all=make.names(group_all)

keep.exprs=filterByExpr(all_data, group=group_all)
data2=all_data[keep.exprs,]

#Voom
x <- calcNormFactors(data2, method = "TMM")
design = model.matrix(~0+group_all)
colnames(design)=make.names(colnames(design))
v=voom(data2, design, plot=TRUE)


#ARMS comparison to Normal Tissue 
contr.matrix = makeContrasts(
  ARMSvsadrenal=group_allARMS - group_allAdrenal.Gland,
  ARMSvsblader= group_allARMS - group_allBladder, 
  ARMSvsblood= group_allARMS - group_allBlood,
  ARMSvsbloodves= group_allARMS - group_allBlood.Vessel,
  ARMSvsbrain= group_allARMS - group_allBrain,
  ARMSvsbreast = group_allARMS - group_allBreast,
  ARMSvscervix= group_allARMS - group_allCervix.Uteri,
  ARMSvscolon= group_allARMS - group_allColon,
  ARMSvsesoph= group_allARMS - group_allEsophagus,
  ARMSvsfallop= group_allARMS - group_allFallopian.Tube,
  ARMSvsheart= group_allARMS - group_allHeart, 
  ARMSvskidney= group_allARMS - group_allKidney,
  ARMSvsliver= group_allARMS - group_allLiver,
  ARMSvslung= group_allARMS - group_allLung, 
  ARMSvsmuscle = group_allARMS - group_allMuscle,
  ARMSvsnerve= group_allARMS - group_allNerve, 
  ARMSvsovary = group_allARMS - group_allOvary,
  ARMSvspancreas = group_allARMS - group_allPancreas,
  ARMSvspituitary = group_allARMS - group_allPituitary, 
  ARMSvsprostate = group_allARMS - group_allProstate, 
  ARMSvssalivary = group_allARMS - group_allSalivary.Gland,
  ARMSvsskin = group_allARMS - group_allSkin, 
  ARMSvssmallint= group_allARMS - group_allSmall.Intestine, 
  ARMSvsspleen = group_allARMS - group_allSpleen, 
  ARMSvsstomach= group_allARMS - group_allStomach, 
  ARMSvstestis = group_allARMS - group_allTestis,
  ARMSvsthyroid = group_allARMS - group_allThyroid, 
  ARMSvsuterus = group_allARMS - group_allUterus, 
  ARMSvsVagina = group_allARMS - group_allVagina,
  levels = colnames(design))


vfit = lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit=eBayes(vfit)
summary(decideTests(efit))


#Sig and fold change 
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


comparison_pvalue=data.frame(row.names = comp3$Gene.name)
comparison_logfc=data.frame(row.names = comp3$Gene.name)
for (i in 1:ncol(tfit)){
  name=colnames(tfit)[i]
  name <- topTreat(tfit, coef=i, n=Inf)
  gene_anot$Transcript.stable.ID= NULL
  gene_anot=unique(gene_anot)
  comp=merge(name, gene_anot, by.x=0, by.y="Gene.stable.ID", all.x=TRUE)
  comp3=comp[order(comp$adj.P.Val),]
  filename=paste(colnames(tfit)[i], ".csv", sep="")
  filepath=paste("/mnt/isilon/diskin_lab/rawan/RNA_seq_R", filename, sep="")
  write.csv(comp3, file = filepath)
  comp_colname=paste(colnames(tfit)[i], "Pval", sep="")
  comp_colname2=paste(colnames(tfit)[i], "Logfc", sep="")
  comparison_pvalue$comp_colname=comp3$P.Value
  colnames(comparison_pvalue)[i]=comp_colname
  comparison_logfc$comp_colname2=comp3$logFC
  colnames(comparison_logfc)[i]=comp_colname2
}

write.csv(comparison_logfc, file="/mnt/isilon/diskin_lab/rawan/RNA_seq_R/comparison_logfc_ARMS.csv")
write.csv(comparison_pvalue, file="/mnt/isilon/diskin_lab/rawan/RNA_seq_R/comparison_pval_ARMS.csv")

