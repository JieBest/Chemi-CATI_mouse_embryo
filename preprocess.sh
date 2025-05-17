
for sample in `cat allsample.txt`; 
do 
  echo $sample
  trim_galore -o $trimout -j 4 --gzip --paired ${inputdir}${sample}_combined_R1.fastq.gz ${inputdir}${sample}_combined_R2.fastq.gz
  rsem-calculate-expression --no-bam-output --paired-end --alignments -p 10 ${rsemout}/${sample}_Aligned.toTranscriptome.out.bam $ref ${rsemout}/${sample}.rsem -q 
done
