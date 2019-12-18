# run_unicycler
Pipeline to run Unicycler (https://github.com/rrwick/Unicycler) assembly on BACs
with  Illumina short reads + PacBio or Oxford Nanopore long reads.

There are several dependency programs required. On the UM Great Lakes Cluster,
assuming you have access to the proper modules, they can be loaded with:

module load picard
module load bowtie2/2.3.4.3
module load racon
module load Pilon
module load spades/3.11.1
module load ncbi-blast/2.9.0
module load samtools
module load minimap2

Other programs are also required to be in your path.  On Great Lakes load the module
KiddLabCustom/2.0.

Required programs include: blat, cutadapt,unicycler.

Java programs should also be accessible with the PICARD_JARS and PILON_JARS environmental variable
set to point to the directory where the .jar files are.  Note, on our cluster this is done by the
module load command.


A test is performed when the program is first executed.

The main program is the driver script run-unicycler.py which does preprocessing
of Illumina reads, filters out hits to E. coli of Illumina and long reads, runs
Unicycler, then rotates the resulting assembly to start in the proper place.  On our 
cluster a sample command would look like:

```
python run_unicycler/run-unicycler.py  \
--name CH82-109H04 \
--outDirBase BAC-assem/CH82-109H04 \
--fqR1 L001_R1_001.fastq.gz \
--fqR2 L001_R2_001.fastq.gz \
--contam eColiK12_DH10B.fa \
--longread fastq_pass/combined-reads.fastq.gz \
--threads 4 --clean --target pTARBAC2_insert.fa 
```
The default value for number of threads is 1, but running with 4 threads seems to work well.
Be sure to modify cluster submission scripts appropriately.

The script rotate-circle.py is also included which can be used to change the start of a
circular assembly.  It can also optionally remove vector sequence after reorienting.

