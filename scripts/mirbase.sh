perl BAM_to_Intermediate.pl -ml 0 -s hsa -pre ../miRbase_Files/hsa.0nt.fasta -gff ../miRbase_Files/hsa.gff2 -mat ../miRbase_Files/mature.fa -bam ../data/mirbase.bam -ref 1 -out test/mirbase.txt -flank 0

perl Intermediate_to_miRspring.pl -in test/mirbase.txt  -s hsa -mat ../miRbase_Files/mature.fa -out test/mirbase.html

