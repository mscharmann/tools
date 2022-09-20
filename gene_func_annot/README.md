# gene functional annotation

after gene structural prediction (e.g. BRAKER2), get a functional annotation

## Inputs
- predicted AA sequences (from e.g. BRAKER2)
- GTF file with genomic coordinates of the predicted coding genes (from e.g. BRAKER2)

## Output
- BLASTP homology searches (outfmt 6) against:
	- Swissprot
	- TAIR
- domain search with InterProscan ( interpro_result.raw )

- aggregate of all the above using a custom python script


## workflow

Paths / filenames are an example..

```
cd /mnt/barn/mscharmann/A_spectabilis_functional_annotation.2022-09

conda create --name func_annot

conda activate func_annot
conda install seqtk snakemake samtools blast gffread bedtools ant -y
conda install -c conda-forge "openjdk=11"
```


get the predicted AAs:
```
cp ../A_spectabilis_BRAKER_annotation.2021-10/wd_braker/augustus.hints.aa ./
```
## Setup BLAST databases
get Arabodopsis genes: (TAIR)
```
wget --no-check-certificate https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_blastsets/Araport_pep_20220119_representative_gene_model.gz
gunzip -c Araport_pep_20220119_representative_gene_model.gz > Araport_pep_20220119_representative_gene_model.pep.fasta
```
the fasta headers have some bad chars; replace them.
```
cat Araport_pep_20220119_representative_gene_model.pep.fasta | sed 's/\xC3/_/g' | sed 's/\xB3/_/g' | sed 's/\x83/_/g' | sed 's/\xC2/_/g' | sed 's/\x82/_/g' | sed 's/\xA0/_/g' | sed 's/\xAD/_/g' | sed 's/\x9F/_/g' > TAIR.fasta
```
make blast database
```
makeblastdb -in TAIR.fasta -dbtype prot -out TAIR
```
cleanup:
```
rm Araport_pep_20220119_representative_gene_model.pep.fasta
```

get SwissProt: https://www.uniprot.org/downloads#uniprotkblink
```
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot
```

## Setup InterProScan
it wants Java 11!
https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html
```
mkdir my_interproscan
cd my_interproscan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.48-83.0/interproscan-5.48-83.0-64-bit.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.48-83.0/interproscan-5.48-83.0-64-bit.tar.gz.md5
```
Recommended checksum to confirm the download was successful:
```
md5sum -c interproscan-5.48-83.0-64-bit.tar.gz.md5
```
Extract the tar ball:
```
tar -pxvzf interproscan-5.48-83.0-*-bit.tar.gz
```
Index hmm models
"Before you run interproscan for the first time, you should run the command: This command will press and index the hmm models to prepare them into a format used by hmmscan."
```
cd interproscan-5.48-83.0/
python3 initial_setup.py
```
Now test running it:
```
./interproscan.sh –-disable-residue-annot –-goterms --applications PANTHER,Pfam --iprlookup --cpu 12 -i test_all_appl.fasta -f tsv -dp
```
InterProscan options to use:
–-disable-precalc -dp  # we do not want to look-up pre-calculated results, because this is a new genome, we are sure that NONE of the sequences we have have ever been submitted before..
–-disable-residue-annot # no need for this level of detail
–-goterms # yes, want them.
--applications PANTHER,Pfam


# Snakefile
```
nchunks = 90
chunks = range(0,nchunks) 
queryfile = "augustus.hints.aa"

rule all:
	input:
		"blastresult.TAIR.txt.gz",
		"blastresult.SProt.txt.gz",
		"interproscan_result.txt.gz"
	shell:
		"""
		# some cleaning up
		if ls chunk_* 1> /dev/null 2>&1
		then
			rm chunk_*
		fi
		if ls fasta.bed* 1> /dev/null 2>&1
		then
			rm fasta.bed
		fi
		"""

rule concat_interproscan:
	input:
		expand( "chunk_{chunk}_interproscan_result.txt", chunk=chunks)
	output:
		"interproscan_result.txt.gz"
	shell:
		"""
		cat {input} | gzip -c > {output}
		"""

rule interproscan:
	input:
		"chunk_{chunk}.fasta"
	output:
		temp("chunk_{chunk}_interproscan_result.txt")
	threads: 2
	shell:
		"""
		## for Interproscan MUST remove "*" stop characters..
		cat {input} | sed "s/\*//g" > {input}.clean
		my_interproscan/interproscan-5.48-83.0/interproscan.sh -dra -dp -iprlookup –goterms -cpu {threads} -f tsv -T tempdir_{input} -i {input}.clean -o {output}
		sleep 2
		rm -r tempdir_{input}
		rm {input}.clean
		"""


rule concat_TAIR:
	input:
		expand( "chunk_{chunk}_blastresult.TAIR.txt", chunk=chunks)
	output:
		"blastresult.TAIR.txt.gz"
	shell:
		"""
		cat {input} | gzip -c > {output}
		"""

rule concat_SProt:
	input:
		expand( "chunk_{chunk}_blastresult.SProt.txt", chunk=chunks)
	output:
		"blastresult.SProt.txt.gz"
	shell:
		"""
		cat {input} | gzip -c > {output}
		"""

rule blast_TAIR:
	input:
		"chunk_{chunk}.fasta"
	output:
		temp("chunk_{chunk}_blastresult.TAIR.txt")
	threads: 4
	shell:
		"""
		blastp -num_threads {threads} -outfmt 6 -evalue 1e-5 -query {input} -db TAIR -out {output}
		"""

rule blast_SProt:
	input:
		"chunk_{chunk}.fasta"
	output:
		temp("chunk_{chunk}_blastresult.SProt.txt")
	threads: 4
	shell:
		"""
		blastp -num_threads {threads} -outfmt 6 -evalue 1e-5 -query {input} -db uniprot_sprot.fasta -out {output}
		"""


rule split_query_file:
	input:
		queryfile
	output:
		expand( "chunk_{chunk}.fasta", chunk=chunks)
	params:
		num_files = nchunks
	shell:
		"""
		samtools faidx {input}
		cat {input}.fai | awk '{{print $1"\t0\t"$2}}' > fasta.bed
		
		# Work out lines per file.
		total_lines=$(wc -l <fasta.bed)
		((lines_per_file = (total_lines + {params.num_files} - 1) / {params.num_files}))
		
		# Split the actual file, maintaining lines.
		split -d --lines=$lines_per_file fasta.bed chunk_
		
		mv chunk_00 chunk_0
		mv chunk_01 chunk_1
		mv chunk_02 chunk_2
		mv chunk_03 chunk_3
		mv chunk_04 chunk_4
		mv chunk_05 chunk_5
		mv chunk_06 chunk_6
		mv chunk_07 chunk_7
		mv chunk_08 chunk_8
		mv chunk_09 chunk_9
		
		# get fasta seqs and remove ">ID:1-XX" appended to SEqID headers by seqtk 
		for x in $(ls chunk_* ); do
			seqtk subseq {input} $x | sed "s/:1-[0-9]*//g" > $x.fasta
		done
	
		"""
```

## Run BLASTs and InterproScan with Snakemake over chunks of the data
```
snakemake --cores 48 --keep-going --rerun-incomplete
```

## Download GAFs, GO terms, other info for the sequences
```
wget http://current.geneontology.org/annotations/tair.gaf.gz

wget http://current.geneontology.org/annotations/goa_uniprot_all.gaf.gz # this is more complete, but takes very long to search; is 16 GB zipped..
wget http://current.geneontology.org/annotations/goa_uniprot_all_noiea.gaf.gz # this is less complete, but VERY quick, only 6 Mbp zipped.


wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go

# PlantTFDB v5.0 Plant Transcription Factor Database
wget http://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz

# link the GTF
ln -s ../A_spectabilis_BRAKER_annotation.2021-10/wd_braker/augustus.hints.gtf ./
```
## Collect all annotations into a TAB file
may have to modify some inut file names in the script ```collect_annotations.py```
```
python collect_annotations.py augustus.hints.aa
```
### final result:
```augustus.hints.aa.functional_annotations.txt```

## clean up
```
rm uniprot_sprot.fasta*  TAIR* interpro2go tair.gaf.gz goa_uniprot*  Araport_pep_20220119_representative_gene_model.gz Ath_TF_list.txt.gz Araport_pep_20220119_representative_gene_model.pep.fasta
rm -rf my_interproscan
```

