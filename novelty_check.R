library(tidyverse)
library(data.table)

## Read in files containing gene boundaries
ens2hgcn <- fread("Ensemble2HGNC.txt.gz")
colnames(ens2hgcn) <- c("ensembl_id", "hgnc")
gene_boundaries <- fread("EnsembleMap.txt.gz", header = T)
colnames(gene_boundaries) <- c("ensembl_id", "ensembl_chromosome", "strand", "ensembl_bp_start", "ensembl_bp_end", "ensembl_id2")

gene_boundaries %>% 
  inner_join(., ens2hgcn, by = "ensembl_id") -> gene_boundaries

## Read novel variants
novel_candidates <- fread("Novel_variants.csv")

window <- 5e5

## Gather all genes in a window around the potentially novel SNP
genes_in_window <- list()
for (i in 1:nrow(novel_candidates))
{
  curr_snp <- novel_candidates[i,]
  
  gene_boundaries %>% 
    filter(ensembl_chromosome == paste0("chr", curr_snp$CHR)) %>% 
    filter(ensembl_bp_start > curr_snp$POS - window & ensembl_bp_end < curr_snp$POS + window) %>% 
    distinct(hgnc, .keep_all = T) %>% 
    mutate(SNP = curr_snp$rsid,
           Original_gene = curr_snp$`Nearest gene`) %>% 
    dplyr::select(Gene=hgnc, SNP, Original_gene) -> genes_in_window[[i]]
  
}

## Save alternative genes to a file for later use
genes_snp_window <- bind_rows(genes_in_window)
genes_snp_window %>% 
  dplyr::select(Gene, SNP) %>% 
  rbind(., novel_candidates %>% dplyr::select(Gene = `Nearest gene`, SNP=rsid)) %>% 
  distinct(Gene) %>% 
  write.table(., "alernative_gene_list.txt",
              col.names = F, row.names = F, quote = F, sep = "")


## Parse scientific articles for mentions of alternative gene names (the nearest gene should also be included)
## The script novelty_check.py produces a file called 'novelty_lookup_genes.results' containing the genes, number of mentions (occurrences) and page number
search_loci <- "python3 novelty_check.py"
system(search_loci)

## Read in the results of the parsing of the scientific articles
novelty_check <- fread("novelty_lookup_genes.results")
novelty_check %>% 
  filter(Occurrences != 0) %>% 
  group_by(Keyword) %>% 
  summarise(Ocurrences = sum(Occurrences),
            Page_no = paste0(`Page Number`, collapse = ",")) %>% 
  arrange(desc(Ocurrences)) %>% 
  inner_join(., genes_snp_window, by = c("Keyword" = "Gene")) %>% 
  inner_join(., novel_candidates, by = c("SNP" = "rsid")) %>% 
  dplyr::select(CHR, POS, SNP, `Nearest gene`, Other_genes=Keyword, Ocurrences, Page_no) %>% 
  arrange(CHR, POS)


## Point to directory containing results files from OpenTargets and GWAS catalogs of nearest genes
ot_gc_files <- "novelty_genes_opentarts_gwascatalog/"
## Search for phenotype of interst (should be case-insenstivie and contain possible alternative uses of the name) if the phenotype is mentioned among locus/gene associations
## I know that no one uses 'fenotype', but this is just an illustrative example
keyword_pheno <- 'phenotype|fenotype'
output_pheno_search <- "novelty_lookup_pheno.results"

## extract lines containing the phenotype of interest
extract_pheno_lines <- paste0("grep -iE 'phenotype|fenotype' ", ot_gc_files, "*.tsv > ", output_pheno_search)
system(extract_pheno_lines)

## Trim the extracted lines
remove_string1 <- "[^/]+_associations_export.tsv:"
remove_string2 <- "[^/]+-associated-studies (1).tsv:"
output_pheno_search_clean <- "novelty_lookup_pheno.results.clean"
trim_file <- paste0("sed -E -e 's|", remove_string1, "||g' -e 's|", remove_string2, "||g' ", output_pheno_search, " > ", output_pheno_search_clean)
system(trim_file)

## Create a gene list for looping through files with phenotype associated SNPs and genes
genes_snp_window %>% 
  dplyr::select(Gene, SNP) %>% 
  rbind(., novel_candidates %>% dplyr::select(Gene = `Nearest gene`, SNP=rsid)) %>% 
  distinct(Gene)  %>% 
  pull(Gene) -> gene_list

## Extract lines containing  the galternative genes and create a file containing 
outpath <- "novelty_check/"
gene_snp_pheno_assoc <- list()
for (gene in gene_list)
{
  # Get the lines containing genes
  get_gene_lines <- paste0("grep ", gene, " ",  output_pheno_search_clean, " | sed -E 's/.*\b(rs[0-9]+[-_a-zA-Z0-9]*)\b.*/\1/g' > ",outpath, gene, "_file.txt")
  system(get_gene_lines)
  
  
  gene_file <- paste0(outpath, gene, "_file.txt")
  
  # read the file if it is non-empty
  if (file.exists(gene_file) && file.info(gene_file)$size > 0) {
    gene_snps <- fread(gene_file, header = F)
    
    gene_snps %>% 
      select(where(~ any(grepl("^rs[0-9]+", .)))) %>% 
      gather() %>%
      mutate(value = gsub("-.*", "", value)) %>% 
      distinct(value) %>% 
      mutate(Gene = gene) %>% 
      dplyr::rename(SNP=value) -> gene_snp_pheno_assoc[[gene]]
    
  }
  
}

## Create a file with 4 columns: phenotype associated SNP, annotated gene, original SNP (novel SNP) and nearest gene to the novel SNP
bind_rows(gene_snp_pheno_assoc) %>% 
  dplyr::rename(Pheno_SNP = SNP) %>% 
  inner_join(., genes_snp_window, by = "Gene") -> snp_compare

## Get LD between phenotype associated SNPs and novel SNPs
## Use populations appropriate for the present study (here we use European population from 1000 genomes)
g1000_freq <- "g1000_eur.afreq"
g1000_bed <- "g1000_eur"
out_ld <- "ld_files/"
curr_snp_compare_r2 <- list()
for (i in 1:nrow(snp_compare))
{
  
  # Get LD between SNP pair
  curr_snp_compare <- snp_compare[i,]
  string1 <- paste0("plink2 --bfile ", g1000_bed, " --ld ", curr_snp_compare$SNP, " ", curr_snp_compare$Asthma_SNP)
  string2 <- paste0(" --read-freq ", g1000_freq, " --out ", out_ld, curr_snp_compare$SNP, "_", curr_snp_compare$Asthma_SNP)
  get_ld <- paste0(string1, string2) 
  system(get_ld)
  
  # Extract r2 from the logfile
  string3 <- paste0("grep r^2 ", out_ld, curr_snp_compare$SNP, "_", curr_snp_compare$Asthma_SNP)
  string4 <- paste0(".log | sed -E 's/[[:space:]]+//g' | sed -E 's/D.*//g' | sed -E 's/^.*=//g' > ", out_ld, curr_snp_compare$SNP, "_", curr_snp_compare$Asthma_SNP, ".ld")
  grep_r2 <- paste0(string3, string4)
  system(grep_r2)
  
  # Insert extracted r2 in the file containing both SNPs
  curr_snp_par <- paste0(out_ld, curr_snp_compare$SNP, "_", curr_snp_compare$Asthma_SNP, ".ld")
  cur_snp_r2 <- fread(curr_snp_par)
  curr_snp_compare %>% 
    mutate(R2 = cur_snp_r2$V1) -> curr_snp_compare_r2[[i]]
  
}

# Save file for comparison (noevl SNPs should be in low LD) with known phenotype-associated SNPs
bind_rows(curr_snp_compare_r2) %>% 
  write.table(., "/home/kasper/nas/closed_projects/polyp_gwas/novelty_check/novelty_lookup_asthma.results.genes.rsid.final",
              col.names = T, row.names = F, quote = F, sep = "\t")



