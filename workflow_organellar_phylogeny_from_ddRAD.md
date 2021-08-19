# workflow for organellar phylogeny from reduced-representation data (RAD, ddRAD, GBS etc)

ddRAD and related methods are now routinely used to build phylogenies, but it is usually not considered that the sequenced loci are actually a mixture of nuclear and organellar loci. As the organellar loci are almost certainly orders of magnitude fewer, their contribution to an indiscriminate multi-locus phylogeny is probably neglegible and the result is essentially a nuclear phylogeny. However, organellar phylogenies can be of special interest. This workflow shows how to use ddRAD (or GBS, or...) reads originating from organellar genomes to generate specific organellar phylogenies. It relies on the availability of an organellar reference assembly for the lineage. Some code snippets are taken from dDocent <https://github.com/jpuritz/dDocent> by J. Puritz. Uses bash, python and R.


dependencies may be installed through conda:

conda install bwa mawk samtools freebayes parallel bedtools bcftools python2.7 raxml r-phangorn r-ape

## 1. read mapping

Prepare a directory with the organelle reference fasta (here named "reference.fasta") and fastq files of the samples (here with prefix "sample_" and suffix ".fq.gz"; some parts of the workflow rely on these pre- and suffices). Then:

```
bwa index reference.fasta

for i in $( ls *.fq.gz | sed 's/.fq.gz//g' ) ; do

bwa mem -a reference.fasta $i.fq.gz -t 6 -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@6 -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
samtools sort -@6 $i.bam > $i
rm $i.bam
mv $i $i-RG.bam
samtools index $i-RG.bam

done

rm *.bam.log bwa.*.log

```

## 2. get mapping success reports
It is good to know whether the reads covered any meaningful portion of the assembly at all.
```
rm mapping_report.txt

all_samples=$(ls ./*.fq.gz | sed 's/.fq.gz//g' )
nbases=$( grep -v ">" reference.fasta | wc | awk '{print $3-$1}' ) 

for sample in $all_samples ; do
	stub_name=$( echo $sample | awk 'BEGIN { FS = "/" } ; { print $2 }' )
	echo $stub_name
	n_raw=$(( $(zcat $sample.fq.gz | wc -l ) / 4 ))	
	n_uniq_mapped=$( samtools flagstat $stub_name-RG.bam | grep "in total (QC-passed reads + QC-failed reads)" | awk 'FS=" +" {print $1} ' )
	prop=$( echo "$n_uniq_mapped / $n_raw" | bc -l )
	n_greater_one=$( samtools mpileup $stub_name-RG.bam | cut -f -4 | awk '{if($4 >= 1) print}' | wc -l )
	n_greater_three=$( samtools mpileup $stub_name-RG.bam | cut -f -4 | awk '{if($4 >= 3) print}' | wc -l )
	n_greater_ten=$( samtools mpileup $stub_name-RG.bam | cut -f -4 | awk '{if($4 >= 10) print}' | wc -l )

	prop_one=$( echo "$n_greater_one / $nbases" | bc -l )
	prop_three=$( echo "$n_greater_three / $nbases" | bc -l )
	prop_ten=$( echo "$n_greater_ten / $nbases" | bc -l )

	echo $stub_name $n_raw $n_uniq_mapped $prop $prop_one $prop_three $prop_ten >> mapping_report.txt
done
```

beautify the mapping reports:
```
cat mapping_report.txt | tr " " "\t" > fu && mv fu ./mapping_report.txt
```
append this header:
sample	raw_reads	mapped_reads	prop_mapped	prop_refbases_cov_greater_0	prop_refbases_cov_greater_2	prop_refbases_cov_greater_9



## 3. call variants
Important: calling HAPLOID genotypes.
```
ls *.bam > inbamlist
freebayes -L inbamlist -v raw.vcf -f reference.fasta --ploidy 1 -m 5 -q 5 -E 3 -G 3 --min-repeat-entropy 1 -V --hwe-priors-off -i -X -u
```


## 4. apply variants to fasta reference & mask missing sites

make a low-coverage .bed file for each sample; here any region with less than 3 mapped reads is detected and will be counted as missing data ("-" character in final alignment). 
 
```
parallel --jobs 20 "genomeCoverageBed -ibam {} -g ./reference.fasta -bga -split | awk '{ if (\$4<3) print}' > {}_regions_below_depth_3.bed" ::: $(ls *.bam)
```
'apply' samples' variants to the reference, and mask missing data (i.e. where insufficient reads were mapped; no data) with "-"
```
reference_genome=reference.fasta
vcf_file=raw.vcf
low_coverage_bed_paths=./

bgzip ${vcf_file}
bcftools index ${vcf_file}.gz

for case in $(bcftools query -l ${vcf_file}.gz ); do
	echo $case
	echo ">"${case} > ${case}_temp.fasta 
	bcftools consensus -m ${low_coverage_bed_paths}/${case}-RG.bam_regions_below_depth_3.bed --fasta-ref ${reference_genome} --missing - -s ${case} ${vcf_file}.gz | grep -v ">" | tr -d "\n" | tr 'N' "-" >> ${case}_temp.fasta
	sed -i -e '$a\' ${case}_temp.fasta
done
```
concatenate samples' fastas to a multi-fasta alignment, cleanup
```
cat sample_*.fasta > ORGANELLE.aln

rm sample_*.fasta
```

## 5. FILTER the alignment
We drop columns from the alignment in which not at least THRESH of all samples have data. Also report missingness per sample.
this is python2.7 code:
```
INFASTA = "ORGANELLE.aln"


indict = {}
with open(INFASTA, "r") as F:
	for line in F:
		if line.startswith(">"):
			seqid = line.strip("\n").strip(">")
		else:
			seq = line.strip("\n")
			indict[seqid] = seq


# how many present sites per sample?
per_sample_missingness_report_dict = {}
for k,v in indict.items():
	present = len([x for x in v if not x == "-"])
	try:
		per_sample_missingness_report_dict[k].append( str((len(v)-present)/float(len(v))) )
	except KeyError:
		per_sample_missingness_report_dict[k] = [ str((len(v)-present)/float(len(v))) ]

		
ntax = len(indict.keys())
thresh = 0.7

min_tax = int(0.7*ntax) 

out_dict = {k:"" for k in indict.keys()}
for col in range(len(indict.values()[0])):
	chars = [x[col] for x in indict.values() if not x[col] == "-"]
	if len(chars) >= min_tax:
		for k in indict.keys():
			out_dict[k] += indict[k][col]

outlines = []
for k in out_dict.keys():
	outlines.append(">"+k)
	outlines.append(out_dict[k])

	
with open(INFASTA+".filtered_min_pres_"+str(thresh)+".fasta.aln", "w") as O:
	O.write("\n".join(outlines)+"\n")

	
outl = ["sample\tmiss_before_filter\tmiss_after_filter"]
for k,v in out_dict.items():
	present = len([x for x in v if not x == "-"])
	per_sample_missingness_report_dict[k].append( str((len(v)-present)/float(len(v))) )
	outl.append(k+"\t"+"\t".join(per_sample_missingness_report_dict[k]))

	
with open(INFASTA + ".per_sample_missingness_report.txt", "w") as O:
	O.write("\n".join(outl)+"\n")


```
How long is the filtered alignment?
```
tail -n 1 *filtered_min_pres_*.aln | wc -c
```

## 6. tree inference
Using RAxML; model GTRGAMMA, 20 ML searches on different random parsimony starting trees, rapid boostrapping 1000. Then make the tree readable for FigTree.

```
rm *.T1 

raxmlHPC-PTHREADS-AVX -T 12 -f a -m GTRGAMMA -p 12345 -x 12345 -# 1000 -s ORGANELLE.aln.filtered_min_pres_0.7.fasta.aln -n T1

cat RAxML_bipartitionsBranchLabels.T1 | perl -p -i -e 's/:(\d+\.\d+)\[(\d+)\]/$2:$1/g' > ORGANELLE.aln.filtered_min_pres_0.7.fasta.aln.GTRGAMMA.BS1000.figtree.tre
```


## 7. collapse short and low-support branches into polytomies 
Explained here: <https://grokbase.com/t/r/r-sig-phylo/138p49bzaz/collapsing-branches-with-low-bootstrap-values>

This R code will export the tree with very short and branches with bootstrap support <70 collapsed into polytomies.

```
library(phangorn)
library(ape)

infilename = "ORGANELLE.aln.filtered_min_pres_0.7.fasta.aln.GTRGAMMA.BS1000.figtree.tre"
outfilename <- paste(infilename,"collapsed_low_support_and_zerobrlen.tre", sep = ".")

intree <- read.tree(infilename)

collapsed_low_support <- pruneTree(intree, 70)

collapsed_zero_brlengths <- di2multi(collapsed_low_support, 0.00001)

write.tree(collapsed_zero_brlengths, file = outfilename, append = FALSE,
           digits = 10, tree.names = FALSE)

```

