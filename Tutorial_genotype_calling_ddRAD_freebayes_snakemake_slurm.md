# Tutorial: genotype calling on ddRAD data with freebayes (Snakemake /w SLURM pipeline)

Freebayes is great but it can be slow and it is not so simple to run it efficiently (= in parallel) on a computer cluster. Snakemake to the rescue!


##### What is Snakemake?
Snakemake is a workflow management system, a tool that helps to automate, organize, and scale complex data analyses. See [Wikipedia](https://de.wikipedia.org/wiki/Snakemake) and  [Github](https://snakemake.github.io/). It has a bit of a steep learning curve... but its really worth it (= saves time and headaches) for any iterative computer cluster stuff. Essentially, you make a list of tasks (possibly interdependent), and Snakemake will execute them for you on the cluster. No need to write any SLURM scripts or so by hand!!

## 0. Inputs
- for each sample, one BAM alignment file with ddRAD data, sorted and indexed (could be created by dDocent mapping step)
- a text file that lists all the BAM files, may create like so: `ls Pop*.bam > bamlist.txt`
- a fasta file with the reference assembly (ddRAD assembly with thousands of small contigs (could be created by dDocent assembly step) 


## 1. Create a conda environment
```
conda create --name gtcall_ddRAD freebayes samtools snakemake -y
```
Then activate the new env like so: `conda activate gtcall_ddRAD`
After this, we still need to install a snakemake plugin, which it requires since version >=8 to interact with cluster job control systems..
(https://stackoverflow.com/questions/77929511/how-to-run-snakemake-8-on-a-slurm-cluster)
```
pip install snakemake-executor-plugin-cluster-generic
```

## 2. Snakefile
Snakemake executes tasks (called "rules") that are written down in a so-called Snakefile. This is a plain text file (use e.g. nano or any other editor), and contains the complete set of instructions of what Snakemake should do. The Syntax of Snakemake is a bit like python, but not quite the same. Create a file `Snakefile` in your working directory and paste the following content into it:

```
import os
import math
import random


reference = "reference.fasta"
bam_list = "bamlist.txt"
freebayes_skip_cov = 20000 # This is a great speed-up and stability parameter for Freebayes. Set this to soemthing like 10 * N_samples * 50 (-> anything that has >500 reads per sample is very likely a collapsed repeat region!)
final_output_name = "test.raw.vcf.gz" # name of the output. MUST end with ".vcf.gz"
chunks = 200 # integer


# Read BAM files
with open(bam_list) as f:
	bam_files = [line.strip() for line in f if line.strip()]

sample_names = [os.path.splitext(os.path.basename(bam))[0] for bam in bam_files]
print(sample_names)
exit
bam_map = dict(zip(sample_names, bam_files))

# Output chunk contig files (flat in root)
contig_chunks = expand("chunk_{i}_contigs.txt", i=range(chunks))

rule all:
	input:
		final_output_name


rule index_ref:
	input:
		reference
	output:
		fai = reference + ".fai"
	shell:
		"""
		samtools faidx {reference}
		"""

rule split_contigs_into_chunks:
	# to clarify, we are NOT splitting any contigs up into smaller contigs, but we make sub-groups of multiple contigs each
	# for better load distribution, contig order (as appears in reference assembly) is randomised before making the chunks. This is because in dDocent, contigs tend to be ordered by the number of reads that were clustered into a contig during the assembly. If left as-is, some freebayes jobs get much more data to crunch than others, and the overall runtime will be longer.
	input:
		fai = reference + ".fai"
	output:
		temp(contig_chunks)
	run:
		contigs = [line.split('\t')[0] for line in open(input.fai)]
		random.shuffle(contigs)
		chunk_size = math.ceil(len(contigs) / chunks)
		for i in range(chunks):
			with open(output[i], 'w') as out:
				out.write('\n'.join(contigs[i*chunk_size : (i+1)*chunk_size]))


rule split_bam_by_contig_group:
	input:
		bam=lambda wc: bam_map[wc.sample],
		contigs=lambda wc: f"chunk_{wc.c}_contigs.txt"
	output:
		bam=temp("{sample}_chunk_{c}.bam"),
		bed=temp("{sample}_bam_chunk_{c}.bed")
	run:
		# Convert contig list to BED file
		with open(input.contigs) as contig_file, open(output.bed, "w") as bed_file:
			for line in contig_file:
				contig = line.strip()
				if contig:
					bed_file.write(f"{contig}\t0\t1000\n")

		shell(
			"samtools view -b -L {output.bed} {input.bam} > {output.bam}"
		)



rule split_reference_by_contig_group:
	input:
		fasta=reference,
		contigs=lambda wc: f"chunk_{wc.c}_contigs.txt"
	output:
		fasta=temp("chunk_{c}.fasta")
	shell:
		"""
		cat {input.contigs} | samtools faidx {input.fasta} -r - > {output.fasta}
		"""



rule call_variants:
	input:
		bams=lambda wc: [f"{sample}_chunk_{wc.c}.bam" for sample in sample_names],
		contigs="chunk_{c}_contigs.txt",
		fasta="chunk_{c}.fasta"
	output:
		temp("chunk_{c}.vcf")
	params:
		bamlist="chunk_{c}_bams.txt"
	run:
		# Write the BAM file list to file
		with open(params.bamlist, 'w') as f:
			for bam in input.bams:
				f.write(f"{bam}\n")

		# ULTRA-important freebayes options for efficiency:
		# -g --skip-coverage N														  
		#		   Skip processing of alignments overlapping positions with coverage >N.
		#		   This filters sites above this coverage, but will also reduce data nearby.
		#		   default: no limit											   
		# https://groups.google.com/g/freebayes/c/5L657BepjpY
		# Erik Garrison writes: "This will jump over sites that have total coverage in all samples >N. This is important when processing genomes with collapsed repeats, 
		# where coverage may go very high. In many genomes, it makes sense to set this to a multiple of your expected read coverage (e.g. 5-10x, or 
		# alternatively several standard deviations above the mean). Variant calls made in these regions of collapsed sequence are unlikely to be "correct" 
		# in that the actual local copy number can be very high, and so the genotype won't be the ploidy that you've set."

		# more explicit of: -m 5 -q 5 -E 3 -G 3 --min-repeat-entropy 1 -V
		#
		#	-m --min-mapping-quality 5 
		#	-q --min-base-quality 10
		#	-E 3 # this is about indels/haplotypes... leave as-is
		#	-G --min-alternate-total 3 # at least this many reads must support an alternate (over all samples) to consider it.
		#	-V --binomial-obs-priors-off   
		shell("""
			freebayes -L {params.bamlist} -v {output} -f {input.fasta} \
			--min-repeat-entropy 1 \
			--min-alternate-total 3 \
			-E 3 \
			--min-mapping-quality 5 \
			--min-base-quality 10 \
			--hwe-priors-off \
			--binomial-obs-priors-off \
			--report-monomorphic \
			--skip-coverage {freebayes_skip_cov}
			
			# now cleanup
			rm {params.bamlist} {input.fasta} {input.fasta}.fai
		""")

rule merge_vcfs:
	input:
		expand("chunk_{c}.vcf", c=range(chunks))
	output:
		final_output_name
	threads: 3
	shell:
		"""
		# MUST NOT USE bcftools concat: it cannot handle non-monotonic POS ordering at interval boundaries.
		# This pipeline is adapted from freebayes-parallel script:
		# https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel

		# Concatenate all chunk VCFs safely
		cat {input} | python $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		"""
```

## 3. set parameters
Section 0 above already lists the BAM files and reference assembly as inputs. These need to be specified in the `Snakefile` itself (this is not the recommended way of doing it but I find it simpler than having an extra YAML file). Modify `Snakefile` on lines 6 and following:


- `reference = "reference.fasta"`
- `bam_list = "bamlist.txt"`
- `freebayes_skip_cov = 20000` # This is a great speed-up and stability parameter for Freebayes. Set this to soemthing like 10 * N_samples * 50 (-> anything that has >500 reads per sample is very likely a collapsed repeat region!)
- `final_output_name = "test.raw.vcf.gz"` # name of the output. MUST end with ".vcf.gz"
- `chunks = 200` # integer value. To let freebayes run in parallel over chunks of the data, the data (i.e. BAM files and reference assembly) will be split in this many Chunks. Then, one freebayes job will be run per chunk, and once all are finished, they will be collected into a single final output VCF.

## 4. RUN!
Running the pipeline is now easy:
- open a Terminal window, ideally with `tmux` so that Snakemake can keep going when you log off from the server. This is fine on a login node, Snakemake does not use much resources.
- essential tmux:
-- Start a tmux session (name will be an integer): `tmux`
-- start named session: `tmux new -s mysession`
-- detach from session: press `Ctrl` + `b` together, then release and press `d`
-- list sessions: `tmux ls`
-- re-attach to session with the name 0: `tmux attach -t 0`
-- close session: re-attach then type `exit`
-- kill all sessions: `tmux kill-server`

- run Snakemake like so, from within the directory where the `Snakefile` is located, and with the conda env active:
```
snakemake -j 100 --restart-times 3 --keep-going --rerun-incomplete --latency-wait 10 \
         --executor cluster-generic \
         --cluster-generic-submit-cmd "sbatch -t 01:00:00 -c 1 --mem=2G"
```
Once this is running, sit back and detach from the tmux session, then check with `squeue` if Snakemake has actually submitted jobs/are they running? Come back after some time to check the progress.

What this does (adjust as needed, **you probably want to increase the job walltime!**):
- `-j100` Snakemake will submit up to 100 jobs or jobs that together use 100 threads in parallel. Whether you actually get 100 parallel jobs depends on SLURM, however.
- `--restart-times 3 --keep-going --rerun-incomplete` When a job fails, it will be re-attempted up to 3 times, the pipeline will continue anyway if possible, and previous output that was not complete will be repeated as well.
- `--latency-wait 10` When an output is missing, Snakemake will wait 10 seconds until it declares the task failed. This is to account for slow disk I/O (writing/reading) or lag in networks
- `--executor cluster-generic` Telling Snakemake to send jobs to a computer cluster instead of running locally on the same machine
- `--cluster-generic-submit-cmd "sbatch -t 01:00:00 -c 1 --mem=2G"` This is the SLURM command that Snakemake will use to send each job to the cluster for scheduling. It here asks for a walltime of 1h, 1 CPU per task, and 2GB of RAM. Adjust as needed, see SLURM syntax e.g. here: https://slurm.schedmd.com/quickstart.html or  https://www.hpc.uni-bayreuth.de/beginners/basic_control/

## 5. DONE!
Check that Snakemake has successfully finished by inspecting the working directory/result VCF.
SLURM jobs will each produce a LOG file, which is mostly not relevant. Cleanup like so: `rm slurm*`



