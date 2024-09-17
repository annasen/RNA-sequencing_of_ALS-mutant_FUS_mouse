# I work in a conda environment with cellranger (cellranger v 7.0.1)
# Cellranger: make fastq

cellranger mkfastq --id=idname_flowcellname \
 --run=/pathway/to/raw/seq/data/ \
 --csv=samples_lane_index.csv --filter-single-index --localcores=15 --localmem=30

# There might be an error: std::exception::what: Barcode lengths in the sample sheet do not match those in --use-bases-mask. That means that the flowcell was sequenced with different index length than what I specified in the sample sheet. To fix this, I could add --use-bases-mask command (check bcl2fastq doc) or add AT in the Sample sheet's barcode, i.e. instead of writing barcode code write the barcode sequence and add AT at the end, 4 lines per sample with 4 different barcodes. List of 10x barcodes: https://github.com/darneson/10XGenomics/blob/master/chromium-shared-sample-indexes-plate.csv

 

# Cellranger: make html report

cellranger count --id=idname \
--fastqs=/pathway/to/fastq/data/ \
--sample=samplename \
--transcriptome=/pathway/to/10x/reference/genome/cellranger/mm10/ --localcores=8 --localmem=64

# If the automatic demultiplexing worked, the respective files are named as the sample name, e.g. Fus.
# If the automatic demultiplexing did not work, there are always two folders, one with fastqs (name of the sample and flowcell) and one with html report (name of the sample and demultiplexed), e.g. wt_flowcellAHMMCCBGXN, wt_demultiplexed



# seqtk

# Because samples were sequenced on different flowcells with different number of cycles, I might need to trim R1 in order to collapse originally different lengths of UMI (e.g. I need 16+12 bp in R1, but got back 26).

seqtk trimfq -e 2 samplename_S1_L001_R1_001.fastq.gz | gzip > samplename_S1_L001_R1_001.fastq.gz
