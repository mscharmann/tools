# build a supermatrix alignment and partition file for genome windows from sequence reads aligned to a reference genome assembly 

Workflow to build a large supermatrix from re-sequecing data (whole-genome libraries, ddRAD, single-digest RAD, GBS etc) mapped to a reference genome; NOT necessarily useful for assemblies made denovo from RAD-like data itself. 
Will use both python2.7 and python3.

further dependencies may be installed through conda:

conda install parallel seqtk bedtools samtools bcftools tabix

for tree building and re-coding add:

conda install raxml r-geiger

## 1. read mapping, variant calling
Use your method of choice for read mapping to the reference, variant calling and filtering. This should result in sorted and indexed .BAM files, and a VCF file with SNPs ONLY. VERY IMPORTANT: This workflow will not work correctly if any INDELs remain in the VCF file! ONLY SNPs can be handled.

Prepare a directory with the organelle reference fasta (here named "reference.fasta"), fastq files of the samples (here with prefix "sample_" and suffix ".fq.gz"; some parts of the workflow rely on these pre- and suffices), and .BAM files with the naming convention "sample_[XXX]-RG.bam". Then:

## 2. create genome windows
```
seqtk comp reference.fasta | awk '{print $1"\t"$2}' > genomefile.txt
bedtools makewindows -w 100000 -g genomefile.txt > windows.bed
rm genomefile.txt
```
## 3. identify missing sites
Some (or many in the case of RAD) sites may not be covered by reads of a sample. These must be excluded to avoid false inferences; if care is not taken at this step, the reference might be filled in for absent data, which is essentially fictitious, synthetic, made-up data...   

Make a .bed file for each sample marking undesired regions; here any region with less than 3 mapped reads is detected and will be counted as missing data ("-" character in final alignment). Adjust number of parallel jobs to your machine. More sophisticated filtering strategies are imaginable; ideally one should filter the invariant sites exactly as filtering the variants (SNPs).. See <https://github.com/ksamuk/pixy> and <https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326>.
 
```
parallel --jobs 20 "genomeCoverageBed -ibam {} -g ./reference.fasta -bga -split | awk '{ if (\$4<3) print}' > {}_regions_below_depth_3.bed" ::: $(ls *.bam)
```
## 4. create fasta alignments per window
For each sample, for each window, 'apply' variants to the reference with bcftools consensus, creating a fasta output. While doing so, mask missing data (i.e. where insufficient reads were mapped; no data) with "-".

Heterozygous variants are exported as IUPAC ambiguity code (bcftools consensus --haplotype I ).

The resulting fasta alignments per window will be placed in a directory "window_fastas".

The raw fasta files can be very large (approx. one reference genome per sample) - to avoid this, immediately discard columns from the fasta alignment which are not present in a minimum number of samples (taxa). We use this small python2.7 script, create and name it "filter_column_presence.py":
```
import sys

thresh = float(sys.argv[1])
INFASTA = sys.argv[2]

#thresh = 0.7
#INFASTA = "window_1.fasta"

inlist = []
with open(INFASTA, "r") as F:
	for line in F:
		if line.startswith(">"):
			seqid = line.strip("\n").strip(">")
		else:
			seq = line.strip("\n")
			inlist.append([seqid, seq])

		
ntax = len(inlist)
min_tax = int(thresh*ntax) 

seqids = [s[0] for s in inlist]

outlist = [ [id,""] for id in seqids ]
for col in range(len(inlist[0][1])):
	chars = [x[col] for x in [s[1] for s in inlist] if not x[col] == "-"]
	if len(chars) >= min_tax:
		for idx in range(len(seqids)):
			outlist[idx][1] += inlist[idx][1][col]

for s in outlist:
	sys.stdout.write(">" + s[0] + "\n" + s[1] + "\n")
	
```

Then run this bash loop below. Adjust the filtering threshold "thresh", it defines which fraction of the total taxa must have non-missing data to retain a column in the alignment. It is usually best to explore trees resulting from complete or near-complete, moderate and very gappy alignments.

```
reference_genome=reference.fasta
vcf_file=FINAL.vcf
low_coverage_bed_paths=./
thresh=0.7

rm *.csi *.vcf.gz

bgzip ${vcf_file}
tabix -f ${vcf_file}.gz

rm -r window_fastas
mkdir window_fastas

while read w ; do
	wname=$(echo $w | awk '{print $1"_-_"$2"_-_"$3}')
	for case in $(bcftools query -l ${vcf_file}.gz ); do
		echo $case
		echo ">"${case} > ${case}_${wname}.temp.fasta 
		echo $w | awk '{print "samtools faidx reference.fasta "$1":"$2+1"-"$3}' > fu ; source fu | bcftools consensus -m ${low_coverage_bed_paths}/${case}-RG.bam_regions_below_depth_3.bed --haplotype I --missing - -s ${case} ${vcf_file}.gz | grep -v ">" | tr -d "\n" | tr 'N' "-" >> ${case}_${wname}.temp.fasta
		sed -i -e '$a\' ${case}_${wname}.temp.fasta
		cat ${case}_${wname}.temp.fasta >> window_${wname}.fasta
		rm ${case}_${wname}.temp.fasta 
	done
	python2.7 filter_column_presence.py ${thresh} window_${wname}.fasta > window_fastas/window_${wname}.filtered.fasta
	rm window_${wname}.fasta
done < windows.bed

```
## 5. concatenate window fastas to a supermatrix

Use the python3 script below:

```
## assemble_supermatrix.py
## python 3

## python assemble_supermatrix.py PATH FILE_SUFFIX MIN_LEN_LOCUS
# reads fasta alignments per-locus from a directory and assemble a supermatrix
# exclude loci that are shorter than MIN_LEN_LOCUS
# name of the fasta file minus FILE_SUFFIX = name of locus in supermatrix model file

import sys, os

allfiles = os.listdir(sys.argv[1])
infiles = sorted([x for x in allfiles if x.endswith(sys.argv[2])])
MIN_LEN_LOCUS = int(sys.argv[3])

print("number of input files:")
print(len(infiles))
print("parsing files to get all taxa")

all_taxa = set()
for f in infiles:
	with open(sys.argv[1] + f, "r") as I:
		for line in I:
			if line.startswith(">"):
				t = line.strip(">").strip()
				all_taxa.add(t)

all_taxa = list(sorted(all_taxa))				
print("taxa:")
print (all_taxa)
print("second pass, now to get sequences and assembling outputs")


locus_count = 0
out_seq_dict = {t:[] for t in all_taxa}
modelfile_lines = []
left_idx = 1
right_idx = 0

for f in infiles:
	nseq = 0
	taxa = []
	seqlen = 0
	with open(sys.argv[1] + f, "r") as I:
		for line in I:
			if line.startswith(">"):
				nseq += 1
				t = line.strip(">").strip()
				taxa.append(t)
			else:
				ll = len(line.strip())
				if ll > seqlen:
					seqlen = ll
#	print(f, nseq, taxa)
	if list(sorted(taxa)) == all_taxa: # taxa are complete
		if nseq == len(taxa): # sequences are complete
			if seqlen >= MIN_LEN_LOCUS:
				print(f)
				locus_count += 1 
				# get the seqs
				with open(sys.argv[1] + f, "r") as I:
					for line in I:
						if line.startswith(">"):
							t = line.strip(">").strip()
						else:
							seq = line.strip("\n")
							out_seq_dict[t].append( seq )
				
				# prep modelfile
				aligned_len = len(seq)
				right_idx += aligned_len
				mline = "DNA" + "," + f.strip(sys.argv[2]) + "=" + str( left_idx ) + "-" + str( right_idx )
				left_idx += aligned_len
				modelfile_lines.append( mline )
				

with open("supermatrix.model","w") as outfile2:
		outfile2.write( "\n".join( modelfile_lines ) + "\n")


for k,v in out_seq_dict.items():
	print(k,len(v),len("".join(v))) 

outlines_fasta = [">" + k + "\n" + "".join(v) for k,v in out_seq_dict.items()]
with open("supermatrix.fasta","w") as outfile1:
		outfile1.write( "\n".join( outlines_fasta ) + "\n")


print("exported " + str(locus_count) + " loci, DONE!")
```

The script will also discard windows (loci) with less than MIN_LEN_LOCUS sites, in the case that windows (loci) will be used to build gene trees. Run the script:

```
python3 assemble_supermatrix.py window_fastas/ .filtered.fasta 100
```

## 6. gene tree inference
Recommend "raxml_for_each_partition_multithreading.py" found here: <https://github.com/mscharmann/tools/tree/master/phylo>. That script will NOT calculate branch supports for the gene trees / partitions.


## 7. Recommended: collapse short branches into polytomies 
RAxML will represent polytomies in NEWICK format as sequential dichotomies with very short branch lengths (1e-6). However, the super-tree tool ASTRAL will only use the tree topology, and will be misled by RAxML's representation of polytomies. Thus, to re-code the polytomies as actual NEWICK polytomies, use "Rscript_collapse_polytomies.txt" found here: <https://github.com/mscharmann/tools/tree/master/phylo>
