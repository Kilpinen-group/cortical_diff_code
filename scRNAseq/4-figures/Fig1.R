## Fig.1b
all <- SetIdent(all, value = all$differentiationDay)
DimPlot(all, label=F, group.by="mid_gest", cols = col_vector[names(col_vector) %in% unique(all$mid_gest)], raster=T) + 
  facet_wrap(all$differentiationDay) + coord_fixed()

## Fig.1f
library(RColorBrewer)
tt <- table(all$differentiationDay, all$mid_unmapped_bucket)
tt <- as.data.frame(round(tt*100/rowSums(tt),2))
colnames(tt) <- c("diffDay","CellType","Percentage")
tt <- as.data.frame(tt)
cols <- col_vector[names(col_vector) %in% tt$CellType]
order <- names(cols)
tt$CellType <- factor(tt$CellType, levels = order)
tt$diffDay <- paste0("Day ", tt$diffDay)

alltp <- sapply(sort(unique(all$differentiationDay)), function(x){
  
  tmp <- all[,all$differentiationDay==x]
  tt <- t(table(tmp$mid_unmapped_bucket, tmp$donorId_simplified))
  tt <- as.data.frame(round(tt*100/rowSums(tt),2))
  tt$tp <- x
  colnames(tt) <- c("Donor","CellType","Percentage","diffDay")
  tt
  
}, simplify=F)

alltp <- do.call("rbind", alltp)
alltp$diffDay <- paste0("Day ", alltp$diffDay, " per donor")
alltp$CellType <- factor(alltp$CellType, levels = order)

tt$Donor <- tt$diffDay
tt$diffDay <- "All donors"
new <- rbind(tt,alltp)

p1 <- ggplot(new, aes(fill=CellType, x=Donor, y=Percentage))+
  facet_wrap(~diffDay, scales = "free_x",nrow=1)+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(plot.title=element_text(size=8, face="bold"),
        axis.text.x=element_text(size=6,angle=90,hjust=1, vjust=0.5),
        legend.title = element_text(size = 8, face="bold", hjust=0.5), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "lines"))+
  guides(shape = guide_legend(override.aes = list(size = 1)),
         fill = guide_legend(override.aes = list(size = 1), ncol=1))+
  scale_fill_manual(name="Cell types",
                    values = cols)+
  scale_y_continuous(labels = scales::percent)+
  xlab("")+
  ylab("Percentage")

p1