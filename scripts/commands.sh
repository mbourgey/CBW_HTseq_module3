mkdir -p $HOME/workspace/HTG/Module3/

docker run --privileged -v /tmp:/tmp --network host -it -w $PWD -v $HOME:$HOME \
--user $UID:$GROUPS -v /etc/group:/etc/group  -v /etc/passwd:/etc/passwd \
-v /etc/fonts/:/etc/fonts/ -v /media:/media c3genomics/genpipes:0.8

export WORK_DIR_M3=$HOME/workspace/HTG/Module3/
export REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/


mkdir -p $WORK_DIR_M3
cd $WORK_DIR_M3
ln -s $HOME/CourseData/HTG_data/Module3/* .

module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/trimmomatic/0.36 mugqic/samtools/1.9 mugqic/bwa/0.7.17 mugqic/GenomeAnalysisTK/4.1.0.0 mugqic/R_Bioconductor/3.5.0_3.7



mkdir -p originalQC/NA12878/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  --threads 2 --regionName ACTL8 --output originalQC/NA12878/


mkdir -p reads/NA12878/

java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12878/NA12878.trim.out


mkdir -p postTrimQC/NA12878/
java -Xmx1G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 2 --regionName ACTL8 --output postTrimQC/NA12878/

mkdir -p alignment/NA12878/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPU:runNA12878_1\tCN:Broad Institute\tPL:ILLUMINA' \
  $REF/genome/bwa_index/Homo_sapiens.GRCh37.fa \
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx2G -jar ${GATK_JAR} SortSam \
  -I /dev/stdin \
  -O alignment/NA12878/NA12878.sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000


#switch to old GATK 3.8
module unload  mugqic/GenomeAnalysisTK/4.1.0.0
module load mugqic/GenomeAnalysisTK/3.8

java -Xmx2G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R $REF/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12878/realign.intervals \
  -I alignment/NA12878/NA12878.sorted.bam \
  -L 1

java -Xmx2G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R $REF/genome/Homo_sapiens.GRCh37.fa \
  -targetIntervals alignment/NA12878/realign.intervals \
  -o alignment/NA12878/NA12878.realigned.sorted.bam \
  -I alignment/NA12878/NA12878.sorted.bam

#return to GATK 4
module unload mugqic/GenomeAnalysisTK/3.8
module load  mugqic/GenomeAnalysisTK/4.1.0.0
java -Xmx2G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  -I alignment/NA12878/NA12878.realigned.sorted.bam \
  -O alignment/NA12878/NA12878.sorted.dup.bam \
  --METRICS_FILE=alignment/NA12878/NA12878.sorted.dup.metrics
#less alignment/NA12878/NA12878.sorted.dup.metrics
java -Xmx2G -jar ${GATK_JAR} BaseRecalibrator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  --known-sites  ${REF}/annotations/Homo_sapiens.GRCh37.dbSNP150.vcf.gz \
  -L 1:17704860-18004860 \
  -O alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -I alignment/NA12878/NA12878.sorted.dup.bam

java -Xmx2G -jar ${GATK_JAR} ApplyBQSR \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -bqsr alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -O alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -I alignment/NA12878/NA12878.sorted.dup.bam

#switch to old GATK 3.8
module unload  mugqic/GenomeAnalysisTK/4.1.0.0
module load mugqic/GenomeAnalysisTK/3.8

java  -Xmx2G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12878/NA12878.sorted.dup.recal.coverage \
  -I alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -L 1:17700000-18100000

#return to GATK 4
module unload mugqic/GenomeAnalysisTK/3.8
module load  mugqic/GenomeAnalysisTK/4.1.0.0

#### Look at the coverage
#less -S alignment/NA12878/NA12878.sorted.dup.recal.coverage.sample_interval_summary
java -Xmx2G -jar ${GATK_JAR} CollectInsertSizeMetrics \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -I alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -O alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv \
  -H alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.histo.pdf \
  --METRIC_ACCUMULATION_LEVEL LIBRARY

#look at the output
#less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv
java -Xmx2G -jar ${GATK_JAR} CollectAlignmentSummaryMetrics \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -I alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -O alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv \
  --METRIC_ACCUMULATION_LEVEL LIBRARY


exit
