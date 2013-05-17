miRspring
=========

Pipeline scripts for creating a miRspring (miRNA sequence profiling) document. 

The product of this pipeline creates a miRspring document, which essentiallyis a small HTML javascript file. 
Each miRspring document has been developed in such a way that it reproduces miRNA sequencing data from a BAM 
file. Additionally it provides unique analysis tools, and  does not need internet connectivity for these processes.


To create a miRspring document the scripts requires:

    - a sorted BAM file and its index, 
    - Samtools.
    - the following miRbase files: hairpin.txt, <species>.gff, mature.fa
    - Perl



To test your system this repository contains small RNA data set (~100,000 reads) for which shell scripts 
have been prepared. Running these shell scripts will create a miRspring document if your system parameters
have been set correctly. 

Please visit the miRspring homepage for detailed information about the project.



[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/davhum/miRspring/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

