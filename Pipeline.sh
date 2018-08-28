ssh -X jd@128.125.247.142

source /home/jd/footprint/bin/activate

htop

################################################################################
# set up/installation
################################################################################
# created conda environment 'wgbs' and installed bismark and bowtie2 (dependencies)
conda create -n wgbs -c bioconda bismark bowtie2 samtools
    # to activate:
        source activate wgbs
          OR
        conda activate wgbs
    #to deactivate
    conda deactivate

# bismark
#   <http://www.bioinformatics.babraham.ac.uk/projects/bismark/>
#   A tool to map bisulfite converted sequence reads and determine
#   cytosine methylation states. The output produced by Bismark discriminates
#   between cytosines in CpG, CHG and CHH context and enables bench scientists
#   to visualize and interpret their methylation data soon after the sequencing
#   run is completed (PMID: 21493656).

# input files supported:
      # sequence format either FastQ or FastA
      # single-end or paired-end reads
      # input files can be uncompressed or gzip-compressed (ending in .gz)
      # variable read length support
      # directional or non-directional BS-Seq libraries
      # full list of alignment modes:
      # <http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_alignment_modes.pdf>

#######################################
# dna-me-pipeline
git clone https://github.com/ENCODE-DCC/dna-me-pipeline.git
    /home/jd/dna-me-pipeline

################################################################################
                        # (1) Index bisulfite genome
################################################################################
# (1) bisulfite convert genome
      # done to allow bowtie alignments
      # need to do only once per genome
# bowtie
      # index bisulfite converted genome
      # refer to bowtie indexing steps

bismark_genome_preparation \
--verbose /media/segil_lab/SegilRaid/staging/Duc/mm10_bis
      # output report:
      # /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/bismark_output.txt
      # reference output:
      # https://hpc.nih.gov/examples/bismark_prep_output.html

################################################################################
                        # download ENCODE data ENCSR803ICQ
################################################################################
cd /media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5

# rep1
wget https://www.encodeproject.org/files/ENCFF810GKB/@@download/ENCFF810GKB.fastq.gz
wget https://www.encodeproject.org/files/ENCFF627KYH/@@download/ENCFF627KYH.fastq.gz
# rep2
wget https://www.encodeproject.org/files/ENCFF110EYU/@@download/ENCFF110EYU.fastq.gz
wget https://www.encodeproject.org/files/ENCFF795IPO/@@download/ENCFF795IPO.fastq.gz

################################################################################
                        # QC fastq files
################################################################################
conda install -c bioconda fastqc
#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
fastqc *.fastq.gz
      # done for rep1 and rep2 fastq.gz files
      output:

conda install -c bioconda trim-galore
      # not sure why this did not work, will use curl instead
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
      # ~/TrimGalore-0.4.5/trim_galore added to alias

trim_galore --paired R1.fastq.gz R2.fastq.gz

# rep1
trim-galore \
--illumina \
ENCFF627KYH.fastq.gz

trim-galore \
--illumina \
ENCFF810GKB.fastq.gz

# rep2
trim-galore \
--illumina \
ENCFF110EYU.fastq.gz

trim-galore \
--illumina \
ENCFF795IPO.fastq.gz


################################################################################
                        # (2) Bismark alignment step
################################################################################
# (2) bismark read alignment
# bisulfite alignment and methylation calling
# input:
      # directory w/ genome of interest
            # folder must have unmodified genoe (.fa)
            # folder must have 2 bisulfite genome subdirectories generated from
            # bismark genome preparation step
      # sequence files to be analyzed (FastQ, FastA)
      # optional parameters
            # note that bismark defaults to directional BS-seq libraries
            # must specify otherwise '--non_directional'
# output:
            # alignment/methylation call output .bam /.sam
            # stats report

# note: talks about strand specific and directionality
# http://www.biostars.org/p/164119/

# sample list:
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5/rep1
      ENCFF627KYH.fastq.gz
      ENCFF810GKB.fastq.gz
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5/rep2
      ENCFF110EYU.fastq.gz
      ENCFF795IPO.fastq.gz

# running bismark
# rep1

bismark \
-u 10000 \
--bowtie2 \
--genome /media/segil_lab/SegilRaid/staging/Duc/mm10_bis \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5/rep1/ENCFF627KYH_trimmed.fq.gz

# running in tmux 11
bismark \
-n 1 \
-l 28 \
--multicore 8 \
--non_directional \
--bowtie2 \
--genome /media/segil_lab/SegilRaid/staging/Duc/mm10_bis \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5/rep1/ENCFF627KYH_trimmed.fq.gz

bismark \
-n 1 \
-l 28 \
--multicore 8 \
--non_directional \
--genome /media/segil_lab/SegilRaid/staging/Duc/mm10_bis \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5/rep1/ENCFF810GKB.fastq.gz

# rep2
# do on backupRaid
sudo rsync -aP /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/Bisulfite_Genome \
/media/segil_lab/BackupRaid/User_Files/Duc

sudo rsync -aP /media/segil_lab/SegilRaid/staging/Duc/mm10_bis/mm10.fa \
/media/segil_lab/BackupRaid/User_Files/Duc

tmux -new-session -s 21 \
#ran this command on backupRaid
# running on tmux 21
bismark \
-u 10000 \
-n 1 \
-l 28 \
--non_directional \
--bowtie2 \
--genome /media/segil_lab/SegilRaid/staging/Duc/mm10_bis \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5/rep1/ENCFF810GKB_trimmed.fq.gz

bismark -u 10000 /REFERENCES/human38/ -1 sample1_F.fq.gz -2 sample1_R.fq.gz

Find ./ -name "ENC*.f*" | xargs rm


# test
bismark \
-u 10000 \
--non_directional \
--bowtie2 \
--genome /media/segil_lab/SegilRaid/staging/Duc/mm10_bis \
/media/segil_lab/SegilRaid/staging/Duc/wgbs/mouse_midbrain_embryo_d10.5/rep2/ENCFF795IPO_trimmed.fq.gz



# tolerates one non-bisulfite mismatch per read
# output files:
#     test_dataset_bismark.bam            --> contains all alignments plus methylation call strings
#     test_dataset_bismark_SE_report.txt  --> contains alignment and methylation summary
# note: current directory must contain sequence files to be analyzed

#option '--pbat'

does not work, try environment with different bowtie2 versions

conda create -n wgbs -c bioconda bismark bowtie2=2.3.0 samtools

conda create -n wgbs_bowtie2.4.0 -c bioconda bismark bowtie2=2.4.0 samtools

################################################################################
                        # post-processing steps
################################################################################
# de-duplication (allows only 1 read for each position in genome)
# to remove PCR artifacts

# filter number of bisulfite conversion related non-bisulfite mismatches

################################################################################
                        # (3) bismark methylation extractor (optional)
################################################################################
# (3) bismark methylation extractor (optional)
      # This step is optional and will extract the methylation information from
      # the Bismark alignment output. Running this additional step allows splitting
      # the methylation information up into the different contexts, getting
      # strand-specific methylation information and offers some filtering options.
      # You can also choose to sort the methylation information into
      # bedGraph/coverage files, or even process them further to genome-wide
      # cytosine methylation reports.

      # can use '--ignore' to ignore bias shown in M-bias plots

# single end
bismark_methylation_extractor -s --gzip sample_bismark_bt2.bam

# paired end
bismark_methylation_extractor -p --gzip sample_bismark_bt2.bam

################################################################################
                        # QC
################################################################################
# M-bias plot
# shows methylation proportoin across each possible position in the read
# can generate graphs with R or excel
# data for M-bias blot written in coverage text file (.cov, .cov.gz)
# format:
<read position> <count methylated> <count unmethylated> <% methylation> <total coverage>

# bismark HTML processing report
bismark2report [options]
# consider doing de-dup

################################################################################
                          # Alignment step
################################################################################
#   map bisulfite sequencing reads against a bismark transformed genome
#
#   dna-me-pipeline v2.1.0
#     <https://github.com/ENCODE-DCC/dna-me-pipeline/releases/tag/v2.1.0>
#     This release implement's parallel processing on the alignment step of the
#     single-end reads pipeline. As run on DNAnexus, this emplys a new
#     dme-align-se-parallel applet which scatters the input reads into a number
#     of fastqs, then aligns each in a separate sub-process, then gathers the
#     resulting alignments into a single bam.
#
#   dme-align-se
#
#   bowtie
#   <https://www.encodeproject.org/software/bowtie/?version=1.0.1>
#
#   samtools
#     <https://sourceforge.net/projects/samtools/files/samtools/>
#     Samtools is a suite of programs for interacting with high-throughput
#     sequencing data. SAMtools implements various utilities for post-processing
#     alignments in the SAM format, such as indexing, variant caller and
#     alignment viewer, and thus provides universal tools for processing
#     read alignments (PMID:19505943).

# mott-trim.py
#     <https://github.com/ENCODE-DCC/dna-me-pipeline>
#     this pipeline then extracts the CpG, CGH, and CHH methylation patterns
#     genome wide. The WGBS pipeline inputs gzipped DNA-sequencing reads (fastqs)
#     and a Bismark-transformed, Bowtie-indexed genome in a tar.gz archive file.
#     These are processed to generate bam alignment files, which in turn produce:
#         Methylation state at CpG (bedMethyl.gz and bigBed)
#         The methylation state at CHG (bedMethyl.gz and bigBed)
#         The methylation state at CHH (bedMethyl.gz and bigBed)
#         Raw signal files of all reads (bigWig)
#         SamTools quality metrics, Bismark quality metrics
#         Pearson correlation, calculated from the two replicates' methylation states at CpG.



# QC raw sequence files
# FastQC
# <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>



sudo rsync -raz --progress /media/segil_lab/SegilRaid/staging/Duc/defcom/ \
/media/segil_lab/Storage/Duc/defcom1
