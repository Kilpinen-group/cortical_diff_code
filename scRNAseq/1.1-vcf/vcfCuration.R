library(GenomicRanges)
library(VariantAnnotation)
library(GenomeInfoDb)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38
bsgenome <- BSgenome.Hsapiens.NCBI.GRCh38


setwd("Genotypes/")

## Load map file
mapFile <- read.table("KIL_JUN_2021.map")

## VCF generated using map and Ped-plus using Plink 1.9
## Load vcf file
vcfFile <- readVcf(file="plink.vcf")

## keep variants in autosomal + X chr
mapFile <- mapFile[!mapFile$V1 %in% c("X","Y","XY","MT","0"),]
vcfFile <- vcfFile[names(rowRanges(vcfFile)) %in% mapFile$V2,]

## Remove unwanted seqlevels
seqlevels(vcfFile) <- seqlevels(vcfFile)[2:23]

## Change plink chromosome 23 to X as standard NCBI notation
#vcfFile <- renameSeqlevels(vcfFile, c("23"="X"))
seqlevelsStyle(vcfFile) <- "NCBI"

## Creation of GRanges from map file
mylocs <- GRanges(Rle(paste0(mapFile$V1)),
                  IRanges(start=mapFile$V4,
                          end=mapFile$V4),
                  strand="*")
seqlevelsStyle(mylocs) <- "NCBI"

#Set names to the GRanges
names(mylocs) <- mapFile$V2

#reannotation of variants using dbSNP version 144
annotated <- snpsByOverlaps(snps, mylocs)

## identify duplicated hits (this indicates the presence of multiallelic variants instead of biallelic ones)
hits <- findOverlaps(mylocs, annotated)

## remove duplicated hits (query)
mylocs2 <- mylocs[-queryHits(hits)[duplicated(queryHits(hits))]]

## reannotation of variants using dbSNP version 144
annotated <- snpsByOverlaps(snps, mylocs2)
hits <- findOverlaps(mylocs2, annotated)

## remove duplicated hits (subject)
if (any(duplicated(subjectHits(hits)))){
  annotated <- annotated[-subjectHits(hits)[duplicated(subjectHits(hits))]]
  hits <- findOverlaps(mylocs2, annotated)
}

#QC: check no duplicates remain
stopifnot(!(any(duplicated(queryHits(hits))) | any(duplicated(subjectHits(hits)))))

#subsample
mylocs2 <- mylocs2[queryHits(hits)]

#retain the non-duplicated and filtered-in variants in the vcfFile
vcfFile <- vcfFile[names(rowRanges(vcfFile)) %in% names(mylocs2),]

##reannotation of variants
newRanges <- rowRanges(vcfFile)
newRanges$RefAnnot <- annotated[subjectHits(hits)]$RefSNP_id


## Need to force seqinfo mix to run snpsByOverlaps
seqinfo(newRanges) <- seqinfo(bsgenome) 
genome(seqinfo(newRanges)) <- "GRCh38.p2"

## Produce the GPos object with the SNPs
annotVar <- snpsByOverlaps(snps, newRanges)

## Infer reference and alternative alleles
refAlt <- inferRefAndAltAlleles(annotVar, bsgenome)
newRanges$inferredREF <- refAlt$ref_allele
newRanges$inferredALT <- refAlt$alt_alleles

## Filter-out variants with multiple alleles
newRanges <- newRanges[elementNROWS(newRanges$inferredALT)==1,]

## About 100K variants of the total 650K do not show the expected reference allele, and are removed.
## That's not a problem as we do not need that many in demuxlet
table(newRanges$REF==newRanges$inferredREF)
# FALSE   TRUE 
# 100808 548661 

newRanges <- newRanges[newRanges$REF==newRanges$inferredREF,]
vcfFile <- vcfFile[names(rowRanges(vcfFile)) %in% names(newRanges),]

## Only 185K variants out of the remaining 550K do show the same REF-ALT than the inferred one.
table(paste0(newRanges$REF,"-",as.character(unlist(newRanges$ALT)))==paste0(newRanges$inferredREF,"-", as.character(newRanges$inferredALT)))
# FALSE   TRUE 
# 363993 184668

## The reason for that is that most of the variants in the vcf are homozygous/unidentified for all the six samples
## and plink cannot inferr the ALT allele
table(as.character(unlist(newRanges$ALT))=="")
# FALSE   TRUE 
# 184668 363993 

table(geno(vcfFile[as.character(unlist(newRanges$ALT))==""])$GT)
# ./.     0/0 
# 1640 2182318

newRanges <- newRanges[as.character(unlist(newRanges$ALT))!="",]
vcfFile <- vcfFile[names(rowRanges(vcfFile)) %in% names(newRanges),]


## Genotype coding allele looks fine (just some missing ones, which demuxlet can handle)
table(geno(vcfFile)$GT)
# ./.    0/0    0/1    1/1 
# 1196 616696 394488  45258


## Let's use reannotated variants ids, instead of Illumina name
names(rowRanges(vcfFile[match(names(newRanges), names(rowRanges(vcfFile))),])) <- newRanges$RefAnnot

#seqlevelsStyle(vcfNew) <- "NCBI"
tmp <- as.data.frame(Seqinfo(genome="GRCh38"))[1:22,]
tmp <- tmp[order(rownames(tmp)),]
seqlevels(vcfFile) <- seqlevels(vcfFile)[order(seqlevels(vcfFile))]

seqinfo(vcfFile) <- Seqinfo(seqnames=rownames(tmp),
                           seqlengths = tmp$seqlengths,
                           isCircular = tmp$isCircular,
                           genome = tmp$genome)

fileData=DataFrame(Value=format(Sys.time(), "%Y%d%m"))
rownames(fileData) <- "fileDate"

fileformat=DataFrame(Value="VCFv4.3")
rownames(fileformat) <- "fileformat"

phasing=DataFrame(Value="none")
rownames(phasing) <- "phasing"

reference=DataFrame(Value="GRCh38")
rownames(reference) <- "reference"

meta(metadata(vcfFile)$header)$source$Value <- "Custom"

meta(metadata(vcfFile)$header)$contig <- DataFrame(length=as.data.frame(seqinfo(vcfFile))$seqlengths,
          row.names=rownames(as.data.frame(seqinfo(vcfFile))))

meta(metadata(vcfFile)$header) <- DataFrameList("fileDate"=fileData,
                                                "fileformat"=fileformat,
                                                "phasing"=phasing,
                                                "reference"=reference)

meta(metadata(vcfFile)$header) <- meta(metadata(vcfFile)$header)[c(3:6)]

## subset only heterozgous variants


vcfName <- "Genotypes/curatedGeno.vcf"
writeVcf(vcfFile,
         file=vcfName,
         index=F)


