# These are some of my scripts that I wrote, with a lot of inspiration from other scripts.

## latest addition: demultiplexing PE RAD data with variable barcode lengths:
readfq_filter_phred20_in_barcode.demultiplex.qual_filter.PE_reads.py

this script replaces process_radtags of the Stacks pipeline
process_radtags has a problem with barcodes of variable length:
- process_radtags can not check barcodes of variable length at the same time
- one may think: solution is to check them one after another, always re-using the discarded reads from previous runs.
	BUT: discarded reads are not output in order => for PE data, pairing information (= read order in the fastq file) is lost!?
	https://groups.google.com/forum/#!msg/stacks-users/37Fq_ucVwRE/4ygo2ocQiUIJ
- suggested solution by J. Catchen: "use process_shortreads, because it preserves order, then afterwards process_radtags to check restriction site"
However, this functionality does not actually seem to work; process_shortreads does NOT preserve the input order of reads, not even the num,ber of fwd and rev reads are the same in the discard files. So I wrote my own solution, which is probably terribly slow, but at least it gets the job done and I know exactly what is going on.