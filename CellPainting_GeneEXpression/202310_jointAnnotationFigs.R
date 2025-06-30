setwd("D:/PostdocUnsync/09_characterizationPaper/202310_jointAnnotationFigs")

library(oob)
library(ggplot2)

source("D:/PostdocUnsync/05_cellPainting/202305_CPfeaturesAnalysis/analyseFeatureMat/function.R")

gexAnnot<-formatAnnotFromMeta(fastRead("inputs/gexAnnot.tsv"),fastRead("inputs/metaAnnotGex.tsv"))
cpAnnot<-formatAnnotFromMeta(fastRead("inputs/cpAnnot.tsv"),fastRead("inputs/metaAnnotCP.tsv"))



colorScales<-c(attr(cpAnnot,"colorScales"),attr(gexAnnot,"colorScales")["finalAnnotation"])

# commonLine<-intersect(unique(gexAnnot$parentalCellLine), levels(cpAnnot$parentalCellLine))
# 
# gexAnnot<-gexAnnot[gexAnnot$parentalCellLine %in% commonLine,]
# cpAnnot<-cpAnnot[cpAnnot$parentalCellLine %in% commonLine,]


gexAnnot$bioRep[is.na(gexAnnot$bioRep)]<-1
gexAnnot$cellLine_replicate<-paste0(gexAnnot$cellLine,"-",gexAnnot$bioRep,".",gexAnnot$techRep)

gexAnnotOfDDay<-gexAnnot[gexAnnot$diffDay=="d20",]
table(gexAnnotOfDDay$sampleId,gexAnnotOfDDay$cellLine_replicate)

parentalCellLinebyLineGex<-sapply(unique(gexAnnot$cellLine), function(line) {
    gexAnnot[gexAnnot$cellLine==line,"parentalCellLine"][1]
} )

parentalCellLinebyLineCP<-sapply(levels(cpAnnot$cellLine), function(line) {
    cpAnnot[cpAnnot$cellLine==line,"parentalCellLine"][1]
} )
names(parentalCellLinebyLineCP)<-levels(cpAnnot$cellLine)

ggGex<-list()
ggCP<-list()
for(dday in c("d20","d40","d70")){
    gexAnnotOfDDay<-gexAnnot[gexAnnot$diffDay==dday,]
    cpAnnotOfDDay<-cpAnnot[cpAnnot$diffDay==dday,]
    
    finalAnnotEnoughCells<-levels(gexAnnotOfDDay$finalAnnotation)[
        table(gexAnnotOfDDay$finalAnnotation) / length(gexAnnotOfDDay$finalAnnotation) * 100 > 1
    ]
    
    gexAnnotOfDDay<-gexAnnotOfDDay[gexAnnotOfDDay$finalAnnotation %in% finalAnnotEnoughCells,]
    
    enrichGex<-enrich2vector.V2(gexAnnotOfDDay$cellLine_replicate,as.character(gexAnnotOfDDay$finalAnnotation), varnames = c("line","cellType"))
    enrichGex$parentalCellLine<-parentalCellLinebyLineGex[enrichGex$line |> as.character()]  |> as.factor()
    
    #ggGex[[dday]]<-
        ggplot(enrichGex,aes(x=log2(OR),y=AMI,fill=parentalCellLine))+
        geom_point(pch=21,size=3)+
        scale_fill_manual(values=colorScales$parentalCellLine[levels(enrichGex$parentalCellLine)])+
        theme_bw()+
        facet_wrap(~cellType)+
        ggtitle(paste0("Day ",dday, " CP"))+
        xlim(c(0,1))
    
    
    CpClusterEnoughCells<-levels(cpAnnotOfDDay$cellClusters)[
        table(cpAnnotOfDDay$cellClusters) / nrow(cpAnnotOfDDay) * 100 > 1
    ]
    
    cpAnnotOfDDay<-cpAnnotOfDDay[cpAnnotOfDDay$cellClusters %in% CpClusterEnoughCells,]
    
    enrichCP<-enrich2vector.V2(cpAnnotOfDDay$cellLine,as.character(cpAnnotOfDDay$cellClusters), varnames = c("line","cellCluster"))
    enrichCP$parentalCellLine<-parentalCellLinebyLineCP[enrichCP$line |> as.character()] |> as.factor()
    
    ggCP[[dday]]<-ggplot(enrichCP,aes(x=log2(OR),y=AMI,fill=parentalCellLine))+
        geom_point(pch=21,size=3)+
        scale_fill_manual(values=colorScales$parentalCellLine[levels(enrichCP$parentalCellLine)])+
        theme_bw()+
        facet_wrap(~cellCluster)+
        ggtitle(paste0("Day ",dday, " CP"))
    
}

pdf("outPlot/abundancePerDDayAMI.pdf",width=20,height=9)
for(dday in c("d20","d40","d70")){
    multiplot(ggGex[[dday]],ggCP[[dday]],cols=2)
}
dev.off()


enrichBoth<-enrichGex
enrichBoth$enrichBoth<-enrichCP$cellCluster


enrichCellType<-enrich2vector(gexAnnot$parentalCellLine,gexAnnot$finalAnnotation)


abundances <- table(gexAnnot$finalAnnot,gexAnnot$cellLine) 
abundances <- unclass(abundances) 

library(edgeR)
# Attaching some column metadata.
extra.info <- gexAnnot[match(colnames(abundances), gexAnnot$cellLine),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$parentalCellLine)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~factor(y.ab$samples$diffDay) + factor(y.ab$samples$parentalCellLine), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

plotBCV(y.ab, cex=1)
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res, n=Inf)



#####AMI plot new ver
gexAnnot<-formatAnnotFromMeta(fastRead("inputs/gexAnnot.tsv"),fastRead("inputs/metaAnnotGex.tsv"))
cpAnnot<-formatAnnotFromMeta(fastRead("inputs/cpAnnot.tsv"),fastRead("inputs/metaAnnotCP.tsv"))

gexAnnot$bioRep[is.na(gexAnnot$bioRep)]<-1
gexAnnot$cellLine_replicate<-paste0(gexAnnot$cellLine,"-",gexAnnot$bioRep,".",gexAnnot$techRep,".",gexAnnot$Protocol)

parentalCellLinebyLineGex<-sapply(unique(gexAnnot$cellLine_replicate), function(line) {
    gexAnnot[gexAnnot$cellLine_replicate==line,"parentalCellLine"][1]
} )

parentalCellLinebyLineCP<-sapply(levels(cpAnnot$cellLine), function(line) {
    cpAnnot[cpAnnot$cellLine==line,"parentalCellLine"][1]
} )


els<-expand.grid(c("Gex","CP"),c("d20","d40","d70")) |> apply(1,function(x) paste0(x,collapse="_"))

res<-vector("list", length(els))
names(res)<-els

for(day in c("d20","d40","d70") ){
    gexAnnotOfDDay<-gexAnnot[gexAnnot$diffDay==day,]
    cpAnnotOfDDay<-cpAnnot[cpAnnot$diffDay==day,]

    finalAnnotEnoughCells<-levels(gexAnnotOfDDay$finalAnnotation)[
        table(gexAnnotOfDDay$finalAnnotation) / length(gexAnnotOfDDay$finalAnnotation) * 100 > 1
    ]
    
    CpClusterEnoughCells<-levels(cpAnnotOfDDay$cellClusters)[
        table(cpAnnotOfDDay$cellClusters) / nrow(cpAnnotOfDDay) * 100 > 1
    ]
    
    gexAnnotOfDDay<-gexAnnotOfDDay[gexAnnotOfDDay$finalAnnotation %in% finalAnnotEnoughCells,]
    cpAnnotOfDDay<-cpAnnotOfDDay[cpAnnotOfDDay$cellClusters %in% CpClusterEnoughCells,]
    
    enrichGex<-enrich2vector.V2(gexAnnotOfDDay$cellLine_replicate,as.character(gexAnnotOfDDay$finalAnnotation), varnames = c("line","cellGroup"))
    enrichGex$parentalCellLine<-parentalCellLinebyLineGex[enrichGex$line |> as.character()]  |> as.factor()
    enrichGex$day<-day
    enrichGex$Modality<-"Gex"
    enrichGex$Protocol <- strsplit(x = as.character(enrichGex$line),split = ".",fixed = T) |> sapply(function(x) x[len(x)])
    
    res[[paste0("Gex_",day)]]<-enrichGex


    
    enrichCP<-enrich2vector.V2(cpAnnotOfDDay$cellLine,as.character(cpAnnotOfDDay$cellClusters), varnames = c("line","cellGroup"))
    enrichCP$parentalCellLine<-parentalCellLinebyLineCP[enrichCP$line |> as.character()] |> as.factor()
    enrichCP$day<-day
    enrichCP$Modality<-"CP"
    enrichCP$Protocol <- "Livesey"
    
    res[[paste0("CP_",day)]]<-enrichCP
}

resDt<-do.call(rbind,res)
resDt$wrapGroup <- paste0(resDt$Modality,"_",resDt$day, "_",resDt$cellGroup)


panel2Plot <- list(
    c(Modality="Gex",cellGroup = "InMGE",day = "d40"),
    c(Modality="Gex",cellGroup = "InMGE",day = "d70"),
    c(Modality="CP",cellGroup = "d20.endoRetNeg",day = "d20"),
    c(Modality="CP",cellGroup = "d40.endoRet",day = "d40")
)

resDtFiltered <- lapply(panel2Plot, function(panel) {
    resDt[resDt$Modality==panel["Modality"] & resDt$cellGroup==panel["cellGroup"] & resDt$day==panel["day"],]
}) 
resDtFiltered<-do.call(rbind,resDtFiltered)

ggplot(resDtFiltered,aes(x=log2(OR),y=AMI,color=parentalCellLine,shape = Protocol))+
    geom_point(size=6)+
    theme_bw()+
    facet_wrap(~wrapGroup)+
    scale_color_manual(values=colorScales$parentalCellLine[levels(resDtFiltered$parentalCellLine)])
ggsave("outPlot/fig5.pdf",width=11,height=9.5)

pdf("outPlot/allGroupDay.pdf",width=20,height=20)
    print(
        ggplot(resDt,aes(x=log2(OR),y=AMI,color=parentalCellLine,shape = Protocol))+
        geom_point(size=6)+
        theme_bw()+
        facet_wrap(~wrapGroup)+
        scale_color_manual(values=colorScales$parentalCellLine[levels(resDtFiltered$parentalCellLine)])
    )
dev.off()
