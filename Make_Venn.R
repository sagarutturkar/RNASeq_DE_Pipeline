library(VennDiagram)

##R script accepts arguments
args <- commandArgs(trailingOnly = TRUE)

list1 = args[1]
list2 = args[2]
name1 = args[3]
name2 = args[4]
output_tag = args[5]
path  = args[6]


image_out = paste(output_tag, "Venn.png", sep ="_")
file_out  = paste(output_tag, "Overlap.tsv", sep ="_")
image_title = paste("Venn diagram -", output_tag, sep =" ")

data1 <- read.table(list1 ,sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")
data2 <- read.table(list2 ,sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")

venn.diagram(list(data1$Gene_ID, data2$Gene_ID), 
             filename = image_out,
             height = 2000,
             width  = 2000,
             resolution = 300,
             imagetype = "png",
             fill = c("#FEAD72", "#FED976"), 
             #alpha=c(0.5,0.5),
             cex = 2.5, cat.cex = 2, cat.pos = 0,
             cat.fontface= 14,
             main = image_title,
             main.just = c(0,-10), main.cex = 2, 
             scaled = FALSE, euler.d = FALSE,
             category.names=c(name1, name2))

common = intersect(data1$Gene_ID, data2$Gene_ID)
DEseq2 = setdiff(data1$Gene_ID, data2$Gene_ID)
EdgeR  = setdiff(data2$Gene_ID, data1$Gene_ID)

output <-  data.frame(Gene_ID = c(common, DEseq2, EdgeR),
                      DE_in_Method = (c(rep("DESeq2 EdgeR", length(common)),
                                      c(rep("DESeq2", length(DEseq2))),
                                      c(rep("EdgeR", length(EdgeR))))))

#Merge DeSeq2 log2FC
output <- merge(output, data1, by = "Gene_ID", all.x = TRUE)
output <- subset(output, select = c(Gene_ID, DE_in_Method, log2FoldChange))
colnames(output)[3] <- paste(name1,"log2FC", sep = "_")

#Merge EdgeR log2FC
output <- merge(output, data2, by = "Gene_ID", all.x = TRUE)
output <- subset(output, select = c(Gene_ID, DE_in_Method, DESeq2_log2FC, logFC))
colnames(output)[4] <- paste(name2,"log2FC", sep = "_")

#Get Average log2FC
output$Average_log2FC = rowMeans(output[,c(3,4)], na.rm = TRUE)

#link Annotations
#Anno = read.table(path, sep = "\t", stringsAsFactors=FALSE, quote = "")
library("readr")
Anno = read_delim(path, delim  = "\t", quote = "")

message("#################################")
message("Dimensions of result table Before Merging are: Row = ", nrow(output), " and Column = ", ncol(output))
message("#################################")

output = merge(output, Anno, by = "Gene_ID", all.x = TRUE)

message("#################################")
message("Dimensions of result table After Merging are: Row = ", nrow(output), " and Column = ", ncol(output))
message("#################################")

#write output to a file
write.table(output, file=file_out, quote=FALSE, sep='\t', row.names = FALSE)


