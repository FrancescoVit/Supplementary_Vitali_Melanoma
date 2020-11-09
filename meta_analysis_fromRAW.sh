## meta analysis with matson

## download data

1) Go to SRA at accession SRP116709; select only the human targeted metagenomic samples and save accession to file 

2) download SRA object
prefetch --option-file SRA_accession_list.txt 

3) extract fastq
fasterq-dump --split-files *.sra

## prepare data

# use submitted reads from Vitali et al. 

# unzip files

gunzip *.fastq.gz

# join paired reads

micca mergepairs -i *_R1.fastq -o ./metaanalysis_assembled_vitali.fastq -l 32 -d 8 -t 6

# use CUTADAPT to select only V4 region

# VITALI: 341F: 5′-CCTACGGGNGGCWGCAG-3′and 805r: 5′-GACTACNVGGGTWTCTAATCC-3′
# MATSON: 515F-806R no sequence reported, using from citation below
# 515F: GTGCCAGCMGCCGCGGTAA
# 806R: GGACTACHVGGGTWTCTAAT
# Caporaso, J., Lauber, C., Walters, W. et al. Ultra-high-throughput microbial community analysis on the Illumina HiSeq and MiSeq platforms. ISME J 6, 1621–1624 (2012). https://doi.org/10.1038/ismej.2012.8

# Trim 806R
/usr/local/bin/cutadapt -g Forward=GTGCCAGCMGCCGCGGTAA --discard-untrimmed -o ./reads_vitali_metaanalysis_F.fastq ./metaanalysis_assembled.fastq >> clipping_forward.txt

# join matson reads

micca mergepairs -i *_1.fastq -p sra_1 -e sra_2 -o ./metaanalysis_assembled_matson.fastq -l 32 -d 8 -t 6

# use CUTADAPT to remove primers

/usr/local/bin/cutadapt -g Forward=GTGCCAGCMGCCGCGGTAA --discard-untrimmed -o ./metaanalysis_assembled_matson_F.fastq ./metaanalysis_assembled_matson.fastq >> filtering_forward_fromMatson.txt
# 136 reads with primer


/usr/local/bin/cutadapt -a Reverse=GGACTACHVGGGTWTCTAAT --discard-untrimmed -o ./metaanalysis_assembled_matson_R.fastq ./metaanalysis_assembled_matson.fastq >> filtering_reverse_fromMatson.txt
# 26 reads with adapter

# NB as the number ar really low, skipping primer trimming on Matson reads

# join both file

cat metaanalysis_assembled_matsen.fastq reads_vitali_metaanalysis_F.fastq > metaanalysis_all.fastq

## OTU picking for meta analysis

# filter N
micca filter -i ./metaanalysis_all.fastq -o ./metaanalysis_all.fasta --maxns 0

# pick OTUs
#micca otu -m denovo_unoise -i ./metaanalysis_all.fasta -o ./denovo_unoise_metaanalysis -t 6 -c

micca otu -m denovo_greedy --id 0.97 -i ./metaanalysis_all.fasta -o ./denovo_greedy_metaanalysis -t 6 -c

# assign taxonomy

export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/

micca classify -m rdp -i ./denovo_greedy_metaanalysis/otus.fasta --rdp-gene 16srrna -o ./denovo_greedy_metaanalysis/taxa.txt




