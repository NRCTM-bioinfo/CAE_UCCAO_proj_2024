#!/bin/bash
#SBATCH -p SVC # partition (queue)
#SBATCH --job-name=CEA
#SBATCH -n 40
#SBATCH --array=1-9
#SBATCH --exclude=node30
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/star.%N.%A_%a.out # STDOUT
#SBATCH -e _log/star.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

id=`sed -n ${SLURM_ARRAY_TASK_ID}p sample.txt`
echo "${id}"

fq_path=.

fq1=${fq_path}/${id}_R1.fq.gz
fq2=${fq_path}/${id}_R2.fq.gz

gtf_file=${path}/GRCm39/gencode/rsem_star_v34/gencode.vM34.annotation.gtf
rsem_path=${path}/GRCm39/gencode/rsem_star_v34/gencode.v34.rsem
star_path=${path}/GRCm39/gencode/star_v34_2.7.9a
out_path=${path}/rnaseq_sham/rsem_v34

echo ${id}
mkdir ${out_path}/${id}
STAR --runThreadN 40 --genomeDir ${star_path} \
  --readFilesCommand zcat \
  --readFilesIn ${fq1} ${fq2} \
  --outFileNamePrefix ${out_path}/${id}/${id}. \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortingThreadN 40 \
  --quantMode TranscriptomeSAM GeneCounts   

rsem-calculate-expression --paired-end --no-bam-output \
  --alignments -p 40 \
  -q ${out_path}/${id}/${id}.Aligned.toTranscriptome.out.bam \
  ${rsem_path} \
  ${out_path}/${id}/${id}.rsem


