Args <- commandArgs(TRUE);
biom_file = Args[1]
tree_file = Args[2]
seq_file = Args[3]
paste("biom_file:",biom_file)
paste("tree_file:",tree_file)
paste("seq_file:",seq_file)
library("phyloseq")
tmp = import_biom(biom_file) 
tree = read_tree(tree_file)

otu.table = as.data.frame(otu_table(tmp))
tax.tab=as.matrix(tax_table(tmp))

id = tax.tab[,7]=="s__"
tax.tab[,7][id]="unclass"
id = tax.tab[,6]=="g__"
tax.tab[,6][id]="unclass"
id = tax.tab[,5]=="f__"
tax.tab[,5][id]="unclass"

nids=paste("OTU",sep = "" ,1:length(rownames(otu.table))) 

id1 = match(rownames(otu.table),rownames(tax.tab))
id2 = match(rownames(otu.table),tree$tip.label)

rownames(otu.table) = nids
rownames(tax.tab) =nids[id1]
tree$tip.label= nids[id2]


ps = phyloseq(otu_table(otu.table,taxa_are_rows = T),tax_table(tax.tab),phy_tree(tree))

seq.table = read.csv(seq_file,header = F,sep = " ")
jishu = seq(1,nrow(seq.table),by=2)
seq.table[jishu,]$V1 ->seq_id
seq_id = as.character(seq_id)
length(unlist(strsplit(seq_id,"[>]")))->leng 
seq_id = unlist(strsplit(seq_id,"[>]")) 
oushu = seq(0,leng,by=2)
rowname = seq_id[oushu]
seq.table[-jishu,] ->seq.table
seq.table$OTU = rowname
rownames(seq.table) =rowname
seq.table = as.data.frame(seq.table[,-2])
otu.table = as.data.frame(otu_table(tmp))
id =match(rownames(otu.table),rownames(seq.table))

seq.table = seq.table[id,]
id =match(rownames(otu.table),rownames(seq.table))
rownames(seq.table) = nids[id]
rownames(otu.table) = nids

tax.tab = as.data.frame(tax.tab)
data.frame = cbind(seq.table,tax.tab,otu.table)
names(data.frame)[1] = "OTU_seq"
names(data.frame)[2] = "OTU_pipeline_ID"
write.csv(data.frame,file = "./data_frame.csv",row.names=T,quote=F)
