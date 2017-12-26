Args <- commandArgs(TRUE);

input_data_dir = Args[1]
output_data_ps = Args[2]
save_data = Args[3]

print(c("数据所在目录是:",input_data_dir))
print(c("拟输出ps文件:",output_data_ps))
library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
suppressMessages(library("dada2"))
library("phyloseq")
suppressMessages(library("DECIPHER"))
suppressMessages(library("phangorn"))

fns <- sort(list.files(input_data_dir,full.names =T))
print(c("以下是需要处理的测序数据"))
l=length(fns)
l=l/2
print(c("Total :",l))
fns
fnFs <- fns[grep("1.fq.gz",fns)]
print(c("以下是paired——file：1.fq.gz"))
fnFs
fnRs <- fns[grep("2.fq.gz",fns)]
print(c("以下是paired——file：2.fq.gz"))
fnRs

file_path<- file.path(input_data_dir,"filtered")
print(c("过滤之后的测序数据是：",file_path))
if(!file_test("-d", file_path)) dir.create(file_path) #建立filter目录

print(c("Now:将paired——end数据分别join在一起"))
filtFs <- file.path(file_path,basename(fnFs)) #filter文件
filtRs <- file.path(file_path,basename(fnRs))
print(c("join Done"))

print(c("Now:Begin trim reads"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=5,
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)
}
print(c("trim Done"))

print(c("Now:将paired——end数据分别join在一起"))
derepFs <- derepFastq(filtFs) #类似merge的过程
derepRs <- derepFastq(filtRs)
print(c("join Done"))

#命名merge 后的文件
sam.names1 <- sapply(strsplit(basename(fnFs),"\\."),`[`,1)
sam.names2 <- sapply(strsplit(basename(fnRs),"\\."),`[`,1)
names(derepFs) <- sam.names1
names(derepRs) <- sam.names2

print(c("Now:利用dada 计算非生物学突变的错误率：beginning dada(多线程) "))
ddF <- dada(derepFs, err=NULL, selfConsist=TRUE,multithread=TRUE)
ddR <- dada(derepRs, err=NULL, selfConsist=TRUE,multithread=TRUE)

print(c("移除非生物学突变：dada again（多线程）"))
dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=TRUE,multithread=TRUE) ##multithread = TRUE,可以加快运行速度。
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=TRUE,multithread=TRUE)
print(c("dada Done"))

print(c("Now:过滤后将测序正链负链merge到一起"))
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
print(c("merge Done"))
###We now merge together the inferred forward and reverse sequences, removing paired sequences that do not perfectly overlap as a final control against residual errors.

####sequence-table,计算the number of times each sequence was observed in each sample
print(c("Now:利用sequence-table计算the number of times each sequence was observed in each sample"))
seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
print(c("sequence-table Done"))

####去除pcr错误引入的嵌合体
print(c("Now:remove pcr过程中产生的嵌合体"))
seqtab <- removeBimeraDenovo(seqtab.all)
print(c("remove Done"))

####assign taxonomy
print(c("Now: assign taxonomy(dadabase is /home/jianglikun/dada2/rdp_train_set_14.fa)"))
ref_fasta <- "/home/jianglikun/dada2/rdp_train_set_14.fa" ###reference
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
print(c("assign Done"))

#####
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

print(c("Now:draw the NJ_tree"))
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
print(c("draw Done"))

ps <- phyloseq(tax_table(as.matrix(taxtab)),otu_table(seqtab,taxa_are_rows = T),phy_tree(fitGTR$tree)) 
print(c("Conggratulations!! ps file is finally generated!!"))
print(c("U can find the ps file in:",input_data_dir/output_data_ps))

save(ps,file=input_data_dir/output_data_ps)
save.image(file=c(input_data_dir,save_data))

print(c("U can find the R_Data file in",c(input_data_dir,save_data)))
