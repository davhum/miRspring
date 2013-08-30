perl BAM_to_Intermediate.pl -ml 0 -s hsa -flank 35 -pre ../miRbase_Files/hsa.35nt.fasta -gff ../miRbase_Files/hsa.gff2 -mat ../miRbase_Files/mature.fa -bam ../data/genome.bam -ref 0 -out test/genome.txt

perl Intermediate_to_miRspring.pl -in test/genome.txt -s hsa -out test/genome.html -flank 35


