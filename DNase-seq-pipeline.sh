# OBJECTIVE: identify differential DHS sites between WT mESC and TKO mESC (DNMT KO's)
### Steps: trim, concatenate fastq's, filter blacklist, align, motif enrichment, differential between WT and TKO

################################################################################
### Variables
################################################################################
pwd="/media/segil_lab/SegilRaid/staging/Duc/dnase-seq/"

wt1_fq="SRR1973528_GSM1657364_DNASE_WT_1_Mus_musculus_DNase-Hypersensitivity.fastq.gz"
wt2_fq="SRR1973529_GSM1657365_DNASE_WT_2_Mus_musculus_DNase-Hypersensitivity.fastq.gz"
tko1_fq="SRR1973530_GSM1657366_DNASE_TKO_1_Mus_musculus_DNase-Hypersensitivity.fastq.gz"
tko2_fq="SRR1973531_GSM1657367_DNASE_TKO_2_Mus_musculus_DNase-Hypersensitivity.fastq.gz"

wt1_trim="${wt1_fq%%.*}_trimmed.fq.gz"
wt2_trim="${wt2_fq%%.*}_trimmed.fq.gz"
tko1_trim="${tko1_fq%%.*}_trimmed.fq.gz"
tko2_trim="${tko2_fq%%.*}_trimmed.fq.gz"

wt_merged="wt_merged.fq.gz"
tko_merged="tko_merged.fq.gz"

wt_sam="DNASE_WT_Mus_musculus_DHS.sam"
tko_sam="DNASE_TKO_Mus_musculus_DHS.sam"

################################################################################
### BWA, Picard
################################################################################
# QC and trim reads
trim-galore \
--fastqc \
$wt1_fq \
$wt2_fq \

trim-galore \
--fastqc \
$tko1_fq \
$tko2_fq \

cat \
$wt1_trim \
$wt2_trim > \
$wt_merged

cat \
$tko1_trim \
$tko2_trim > \
$tko_merged

# building genome index for alignment
# done once only
bowtie-build \
/media/segil_lab/SegilRaid/Data_Processing/mm9/mm9/mm9.fa \
/media/segil_lab/SegilRaid/Data_Processing/mm9/mm9/mm9.index

# fastq concatenate and alignment
bowtie -p4 --verbose -v 3 -m 1 --best --strata -S \
/media/segil_lab/SegilRaid/Data_Processing/mm9/mm9/mm9.index \
"${wt_merged}" \
"${wt_sam}"

bowtie -p4 --verbose -v 3 -m 1 --best --strata -S \
/media/segil_lab/SegilRaid/Data_Processing/mm9/mm9/mm9.index \
"${tko_merged}" \
"${tko_sam}"

################################################################################
### Hotspot, Samtools - ID regions with signal
################################################################################
samtools view -S -b DNASE_WT_Mus_musculus_DHS.sam > DNASE_WT_Mus_musculus_DHS.bam
samtools view -S -b DNASE_TKO_Mus_musculus_DHS.sam > DNASE_TKO_Mus_musculus_DHS.bam

bedtools bamtobed -i DNASE_WT_Mus_musculus_DHS.bam > DNASE_WT_Mus_musculus_DHS.bed
bedtools bamtobed -i DNASE_TKO_Mus_musculus_DHS.bam > DNASE_TKO_Mus_musculus_DHS.bed

# peak calling
fseq -f 0 -l 51 -of bed -v DNASE_WT_Mus_musculus_DHS.bed
fseq -f 0 -l 51 -of bed -v DNASE_TKO_Mus_musculus_DHS.bed

# depth normalization

ma plot - log2 of mean expression of two conditions
deseq

# signal generation


# file format conversion


################################################################################
### differential signal enrichment between conditions
################################################################################






################################################################################
### motif analysis of TF enrichment
################################################################################







################################################################################
### Additional Analysis
################################################################################
try DNaseR
