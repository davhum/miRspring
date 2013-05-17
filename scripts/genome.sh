perl BAM_to_Intermediate.pl -ml 0 -s hsa -pre ../miRbase_Files/hsa.35nt.fasta -gff ../miRbase_Files/hsa.gff2 -bam ../data/genome.bam -ref 0 -out test/genome.txt
perl Intermediate_to_miRspring.pl -in test/genome.txt -s hsa -gff ../miRbase_Files/hsa.gff2 -mat ../miRbase_Files/mature.fa -pre ../miRbase_Files/hsa.35nt.fasta -out test/genome.html

