# Tutorial: Filter a VCF from a ddRAD-seq dataset with denovo assembly

Raw sequencing data and raw VCF genotype calls are FULL of errors and artefacts, especially denovo-assembled ddRAD-seq datasets. Errors and artefacts are problematic because they are often not immediately apparent and can be positively misleading. But one can also ruin a dataset by throwing away too much good information. Therefore, careful filtering is required.

**The aim is to get a dataset that has as few errors and artefacts as possible, while retaining as many high-quality genotypes (non-empty entries in the VCF) as possible. We want to retain bi-allelic variants AND monomorphic/invariant sites.** The latter are essential for estimates of genetic diversity, estimates of effective population sizes and to scale demographic inference etc.
The VCF must also not be too empty, because many population genetic analyses will give misleading, biased results when information is missing and not controlled for.

The missingness of the dataset is simply: number of missing genotypes divided by number of potentially present genotypes (there is some helper code below to count this in a VCF)

##### What is a VCF file?
VCF stands for Variant Call Format, see [Wikipedia](https://en.wikipedia.org/wiki/Variant_Call_Format) and exact specifications [here](https://samtools.github.io/hts-specs/) (but note that many tools using VCFs are not actually completely following these specifications).


##### example VCF
Here, it is assumed that the VCF was produced after assembling and mapping reads with **dDocent** (Puritz et al. 2014), then ran **freebayes**.

create a list of the input BAM files
````
ls Pop*.bam > list_of_bamfiles.txt
````
run freebayes, for speed we skip sites which over all 14 samples in this example, summed up have more than 3000 reads (3000/14 = 214 reads per sample)

   -g --skip-coverage N
                   Skip processing of alignments overlapping positions with coverage >N.
                   This filters sites above this coverage, but will also reduce data nearby.
                   default: no limit
````
freebayes -L list_of_bamfiles.txt -v raw.vcf -f reference.fasta -m 5 -q 5 -E 3 -G 3 --min-repeat-entropy 1 -V --hwe-priors-off -i -X -u --report-monomorphic --skip-coverage 3000
````

Tutorial is inspired by: 
Puritz JB, Hollenbeck CM, Gold JR. 2014. dDocent: A RADseq, variant-calling pipeline designed for population genomics of non-model organisms. PeerJ 2:e431. doi:10.7717/peerj.431
[dDocent website](https://ddocent.com//)


## About the ideas behind the filtering steps

#### Quality filters

##### genotype filters (single fields in the 2-dimensional VCF matrix)
- min_read_depth, applying same threshold globally. Lower read depth often means a higher probability of errors / artefacts in the genotype.
- max_read_depth, threshold calculcated specifically for each individual sample. 
	Excessively high read depth is a hallmark of artefacts (i.e. repetitive sequences mapping to the same reference genome coordinates ("collapsing")).
	Because samples often vary in the mean read depth simply due to the number of raw reads that were generated (this cannot be precisely controlled due to technical reasons, e.g. wetlab), 
		it is a good idea to have a specific max_read_depth threshold for each sample.
		With a simple global threshold, too much crappy data is retained in samples with relatively low overall depth, and too much good data is removed in samples with relatively high overall depth.

##### site filters (rows in the VCF)
- allele balance. Allele balance is the ratio of the count of observations for the two alleles in heterozygotes. 
	In real heterozygotes, both alleles should be +- equally common in the read data, allele balance ~ 0.5, but we accept 0.25 - 0.75 to give some room for randomness.). Stronger deviation from this expectation is a hallmark of artefacts (i.e. repetitive sequences mapping to the same reference genome coordinates ("collapsing")). 
- heterozygote_excess. Heterozygote excess is often a hallmark of artefacts (i.e. repetitive sequences mapping to the same reference genome coordinates ("collapsing")). BUT this filter must not be used when you are looking for regions that are expected to be heterozygous in many samples, e.g. males with XY sex chromosomes!

#### Missingness filters

We want to maximise both the number of sites in the dataset and the number of individuals, but there is a trade-off.
- To keep all sites (VCF rows) is only possible by throwing out any samples with missing genotypes; these can be most or even all samples.
- To keep all samples (VCF columns) is only possible by throwing out any sites with missing genotyopes; these can be most or even all sites.
 
The compromise is to throw out only the sites with unusually high missing genotypes, and those samples with unusually high missing genotypes. The end result should be a dataset that retains most sites and most samples, while the global proportion of missing data is not zero but acceptable. **Because one often does not know a priori what is "acceptable", it is good practice to produce several datasets with different missingness filtering, and explore how they perform in downstream analyses.** I recommend to explore both the high and low extremes of missingness and one intermediate dataset.


##### site missingness filter (rows in the VCF)
- missingness (proportion of absent genotypes / no information, per site). Throws away sites with too little data.

##### individual missingness filter (columns in the VCF)
- missingness. In every dataset, there are usually some samples that just do not have enough raw data (this cannot be precisely controlled due to technical reasons, e.g. wetlab).



## Filtering steps

We start with a VCF named `raw.vcf`. To run the below, you will need `bcftools`, `vcftools`, `vcflib`, and `python3` (may install with conda).

##### 1. decompose complex variant calls into phased SNP and INDEL genotypes
and keep the INFO flags for loci and genotypes, then remove phased "|" annotation to unphased "/"
```
vcfallelicprimitives raw.vcf --keep-info --keep-geno | sed 's/|/\//g' > step1.vcf
```
##### 2. remove indels, and sites with more than 2 alleles (keep only invariant or bi-allelic sites)
```
vcftools --vcf step1.vcf --remove-indels --max-alleles 2 --recode --recode-INFO-all --stdout > step2.vcf
```

##### 3. maximum depth filter

First write geno-depth output (report of the depths at each genotype)
```
vcftools --vcf step2.vcf --geno-depth
```

Then calculate upper thresholds and filter. Download the helper python script here: `wget https://raw.githubusercontent.com/mscharmann/tools/refs/heads/master/vcf_drop_genotypes_exceeding_3std_indiv_mean_depth.v3.py`
```
cat step2.vcf | python vcf_drop_genotypes_exceeding_3std_indiv_mean_depth.v3.py --gdepth out.gdepth > step3.vcf
```		
##### 4. minimum depth filter
```
vcftools --vcf step3.vcf --minDP 3 --recode --recode-INFO-all --stdout > step4.vcf
```


##### 5. allele balance filter
Find and mark the bad ones with vcffilter, create a BED-format file; structure `contig	0	1000` lets bcftools exclude the contig completely (because they are all shorter than 1000 bp) 
```
vcffilter -F AB_filter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" step4.vcf | grep AB_filter | tail -n +2 | cut -f1 | sort | uniq | awk '{print $1"\t0\t1000"}' > contigs_failing_allele_balance.bed
```
now remove them all
```
bcftools view -T ^contigs_failing_allele_balance.bed -Ov step4.vcf > step5.vcf
```

##### 6. throw away contigs that have any sites with excessive heterozygosity
First create a BED-format file; structure `contig	0	1000` lets bcftools exclude the contig completely (because they are all shorter than 1000 bp)
Hardy-Weinberg test implemented in VCFtools, yields a p-value for excess of heterozygotes in column 8. We find those variants with p<=0.05, and extract a list of the contigs with these offending variants. Then throw all sites on these contigs away (don't want to trust anything on that contig if there is any offending variant on it).
```
vcftools --vcf step5.vcf --hardy
cat out.hwe | awk '{if($8<=0.05) print $1}' | sort | uniq | awk '{print $1"\t0\t1000"}'> contigs_with_heterozygote_excess.bed
bcftools view -T ^contigs_with_heterozygote_excess.bed -Ov step5.vcf > step6.vcf
```


##### 7. remove individuals with high missingness
First create the report and inspect it
```
vcftools --vcf step6.vcf --missing-indv
cat out.imiss
```

Now make a list of the bad samples, `high_missingness_individuals.txt`, and exclude them. May use a text editor like `nano`.
```
vcftools --vcf step6.vcf --remove high_missingness_individuals.txt --recode --recode-INFO-all --stdout > step7.vcf
```

##### 8. cleanup of removed ALT alleles
just book-keeping: and clean up ALT alleles that still appear in the ALT vcf column but are not actually occurring in the genotype fields

ugly workaround for a bcftools 'bug' as hinted towards here: `https://www.biostars.org/p/9523322/`
```
bcftools view -h step7.vcf | sed 's/Number=R/Number=./' | sed 's/Number=A/Number=./' | sed 's/Number=G/Number=./' > fixed_header.txt
bcftools view -H step7.vcf > body.txt
cat fixed_header.txt body.txt > fixed.vcf
rm fixed_header.txt body.txt
```
Now finally do that...
```
bcftools view --trim-alt-alleles fixed.vcf > step8.vcf
rm fixed.vcf
```
##### 9. remove sites with excessive missingness
**a. extremely lax:** keep sites that have non-missing genotypes in at least 50% of the samples
```
vcftools --vcf step8.vcf --max-missing 0.5 --recode --recode-INFO-all --stdout > final_filtered.maxmiss_0.5.2025-07-24.vcf
```

**b. intermediate:** keep sites that have non-missing genotypes in at least 90% of the samples
```
vcftools --vcf step8.vcf --max-missing 0.9 --recode --recode-INFO-all --stdout > final_filtered.maxmiss_0.9.2025-07-24.vcf
```

**b. extremely strict:** keep sites that have no missingness at all, present in 100% of the samples
```
vcftools --vcf step8.vcf --max-missing 1.0 --recode --recode-INFO-all --stdout > final_filtered.maxmiss_1.0.2025-07-24.vcf
```

### DONE!

**Question at the end: How much of the raw data has made it to the end, which steps remove most data?**


## Useful helpers

#### To count the missingness at genotype level from a VCF. This can be pasted into a terminal:
```
vcf=XXX.vcf
bcftools query -f '[%GT\t]\n' $vcf | \
tr '\t' '\n' | \
awk '
BEGIN { total=0; missing=0 }
{
    total++
    if ($1 == "./." || $1 == ".|." || $1 == "." ) {
        missing++
    }
}
END {
    print "Total potential genotypes:", total
    print "Missing genotypes:", missing
    print "Proportion missing:", missing/total
}'
```

#### To count the number of SNPs retained (--mac 1 => minor allele count 1, i.e. invariants fail this filter):
```
vcftools --vcf XXX.vcf --mac 1 
```


#### To visualize the first N sites from a vcf with R. 
- This needs a few R packages, can also be installed with conda like so: `conda install r-ggplot2 r-reshape2 r-vcfr`.
- Create (maybe using `nano`) a file named `plot_genotypes.R` with the content below (chatgpt helped to make it).
- Run it like so: `Rscript plot_genotypes.R raw.vcf 10000`

```
# ---- Load Required Packages ----
suppressMessages(library(vcfR))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(tools))  # for file_path_sans_ext()

# ---- Handle Command Line Arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_genotypes.R <vcf_file.vcf[.gz]> [num_sites]", call. = FALSE)
}

vcf_file <- args[1]
max_sites <- ifelse(length(args) >= 2, as.numeric(args[2]), 10000)
if (is.na(max_sites) || max_sites < 1) stop("Second argument must be a positive integer.")

# ---- Derive Output Filename ----
base <- basename(vcf_file)
base <- sub("\\.vcf(\\.gz)?$", "", base)
output_png <- paste0(base, "_", max_sites, "_genotype_heatmap.png")

# ---- Read VCF ----
cat("Reading VCF:", vcf_file, "\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# ---- Filter to Biallelic Sites (Keep Monomorphic Too) ----
biallelic_sites <- !grepl(",", vcf@fix[, "ALT"])  # exclude multi-allelics
vcf <- vcf[biallelic_sites, ]

# ---- Limit to First N Sites ----
vcf <- vcf[1:min(nrow(vcf), max_sites), ]

# ---- Extract Genotype Matrix ----
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# ---- Convert GTs to Numeric ----
convert_gt <- function(x) {
  ifelse(x %in% c("0/0", "0|0"), 0,
  ifelse(x %in% c("0/1", "1/0", "0|1", "1|0"), 1,
  ifelse(x %in% c("1/1", "1|1"), 2, NA)))
}
gt_numeric <- apply(gt, c(1, 2), convert_gt)

# ---- Reshape for ggplot ----
gt_df <- melt(gt_numeric)
colnames(gt_df) <- c("Site", "Sample", "Genotype")
gt_df$Genotype <- as.factor(gt_df$Genotype)

# ---- Plot ----
p <- ggplot(gt_df, aes(x = Sample, y = Site, fill = Genotype)) +
  geom_tile() +
  scale_fill_manual(
    values = c("0" = "blue", "1" = "purple", "2" = "red"),
    na.value = "grey80",
    name = "Genotype",
    labels = c("HomRef", "Het", "HomAlt")
  ) +
  labs(
    title = paste("Genotype Heatmap:", base, "-", nrow(gt_numeric), " sites"),
    x = "Sample",
    y = "Variant Site"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ---- Save PNG ----
ggsave(
  filename = output_png,
  plot = p,
  width = 8.27,
  height = 11.69,
  dpi = 300,
  units = "in"
)

cat("Genotype heatmap saved to", output_png, "\n")

```



