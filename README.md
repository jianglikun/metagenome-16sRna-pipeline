# 16S_to_biom_pipeline

This script was write to analysis 16s sequence data. Finaly get the sample_group biom file.

Use software: Qiime

The main script:
    1.join_paired_ends.py:    This script takes forward and reverse Illumina reads and joins them using the method chosen

    2.split_libraries_fastq.py:   qulity control for the sequence

    3.pick_open_reference_otus.py:   This script will produce an OTU mapping file (pick_otus.py), a representative set of sequences (FASTA file from pick_rep_set.py), a sequence alignment file (FASTA file from align_seqs.py), taxonomy assignment file (from assign_taxonomy.py), a filtered sequence alignment (from filter_alignment.py), a phylogenetic tree (Newick file from make_phylogeny.py) and a biom-formatted OTU table (from make_otu_table.py).This script will produce an OTU mapping file (pick_otus.py), a representative set of sequences (FASTA file from pick_rep_set.py), a sequence alignment file (FASTA file from align_seqs.py), taxonomy assignment file (from assign_taxonomy.py), a filtered sequence alignment (from filter_alignment.py), a phylogenetic tree (Newick file from make_phylogeny.py) and a biom-formatted OTU table (from make_otu_table.py).

Input: data_path,in this path include the paired_end sequences file.  
        
    eg:1.1.fq.gz; 1.2.fq.gz   they must End with "***.fq.gz"

Otuput:
    Biom file.  
    The biom file is outputfile/otus/otu_table_mc2_w_tax_no_pynast_failures.biom.
    The otu sequence is outputfile/otus/rep_set.fna
    The phyloseq tree is outputfile/otus/rep_set.tre
    
    DATA_table:
    The data_table is data_frame.csv.

other :
    Join paired-end data method:SeqPrep - (https://github.com/jstjohn/SeqPrep)
    
    In pick_open_reference_otus.py: otu_picking_method is usearch61(http://www.drive5.com/usearch/download.html)
    
    In pick_open_reference_otus.py need a parameter file : in this script the parameter file is /home/jianglikun/16s_pipeline/p_file
    
    In the step five i use a R_script /home/jianglikun/16s_pipeline/R_for_biom.R to get the data_table.
