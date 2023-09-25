##### Hi there 👋 😁

I'm a PhD student in Cameroon.
I'm currently working on evolutionary genomic to identify signatures of positive selection associated to resistance escalation in major malaria vectors.
Shell and R scripts below described the pipeline used for the analyses and visualization of the Poolseq GWAS data.

# Quality control script using fastqc and multiqc
I am going to show here how to use the `Fastq_Quality_check.sh` script i wrote to sequentialy perform quality control of your NGS data using `fastqc` then pipe the sdout from `fastqc` into `multiqc` to aggregate the results and visualize. This command helps you to save more time and speed your analyses.
Before to start please ensure that you make the script executable using the following command `chmod +x Fastq_Quality_check.sh`.

To to that,
- invoke `Fastq_Quality_check.sh` then run this command:
  
`./Fast_Quality_check.sh -i path/to/the/folder/containing/FastQ_files -o path/to/the/output_file -t 20`

Where,
- -i represents the path to the folder where all your fastQ files are located;
- -o represents the path to the file where you want to save your results (note that if it doesn't exist, it will be automatically created);
- -t is the number of thread you want to run your job with (you can adjust according to the performance and parameters of your computer).
  
# Alignment script
The `Alignment.sh` script was written to help align samples that have 2 or more read pairs (two or more reads for forward and for reverse reads post pair-end WGS). The alignment uses `bwa mem` algorith to perform alignment and convert the sam output into bam output using samtools.
In details,
- It check the number of forward and reverse reads for each sample.
- If a sample has exactly one pair of forward and reverse reads, the standart aligment usin BWA-MEM will be run;
- If a sample has two, three, or four pairs of forward and reverse reads, it will combine all the forward reads and all the reverse reads into two separate files and then run the BWA-MEM alignment on the combined reads;
- After the alignment is complete, the combined read files will be clean;
- If a sample doesn't have the correct number of reads, a message indicating that the sample is skipped will be print.
- 
How to run the command?
- invoke `Alignment.sh`

`./alignment.sh -b /path/to/base_dir -t 60 -r /path/to/reference/genome.fasta -o /path/to/output_dir`

Where
- -b represents the base directory;
- -t is the number of threads to use for the alignment;
- -r refers to the reference genome, note that it must be indexed before running (`bwa index genome.fasta`);
- -o refers to the sorted bam output file.

# Extraction of unmapped reads for Metagenomic analyses (EUMA)

I am going to show you how to perform metagenomic analyses using unmapped reads from WGS with the aim to identify microbial signature that could be potentially associated to insecticide super-resistance in malaria vectors.
During the WGS of malaria vectors, beside the entire genome that will be sequenced, other organisms genomes present in the sample are also sequenced. This involved bacteria, parasites, Eukaryotes, Viruses, funga and even unknown organisms refer as `dark reads` that will need further investigation to discover new species using De Novo assembly for example.
We are going to use here the script `Extract_unmapped_reads.sh` to retrieve unmapped reads from bam files then convert them into fastq files.
How to run the command?
- invoke `Extract_unmapped_reads.sh` and run

`./Extract_unmaped.sh -i /path/to/input_directory -o /path/to/input_directory/Unmapped_reads -t 20`

Where
- -i represents the input directory containing all the bam files;
- -o represents the output directory where to store all the unmapped reads for each sample;
- -t represents the number of threads for parallel processing (it depend on your computer capacity).
  

# Sorting and marking duplicates from .bam files
Here, i will show you how to sort according to coordinates, mark and remove duplicates using picard tools with bam files as input.
This session will use the shell script `Sorting_marking_duplicates.sh`.

How to run the command line?

`./Sorting_marking_duplicates.sh -i /path/to/input_dir -o /path/to/output_dir -m metrics.txt -p /path/to/picard.jar`

Where,
- -i represents the file in which the bam files to process are into;
- -o represents the file in which the output will be redirected;
- -m represents the file in which the metric files of each bam file will be save after marking duplicates;
- -p represents the path to the tool you will use to process the sorting and deduplication (in this case, picard.jar from the Picard toold).

# Mapping and coverage statistics
Here, i will show you how to quickly compute mapping and coverage statistics using picard tools and samtools. This part will use `Mapping_statistics.sh` and `Coverage_statistics.sh` shell scripts.
How to run the command line for mapping statistics?

`./Mapping_statistics.sh -b path/to/the/bam_files/directory -o path/to/the/output_directory -p /path/to/the/picard/Software/directory/ -r /path/to/the/reference/VectorBase-61_AfunestusFUMOZ_Genome.fasta`

Where,
- -b represents the file in which the aligned bam files are located;
- -o represents the file in which the output will be redirected (it will be created if it doesn't exist);
- -p represents the directory in which the picard.jar is located;
- -r represents the path to the reference genome of the organism of interest.

How to run the command line for coverage statistics ?

`./Coverage_statistics.sh -b /path/to/bam_directory -o /path/to/output_directory`

Where,
-  -b represents the bam files directory;
-  -o represents the output directory (it will be automatically created if it doesn't exist).


# Conversion of fractions into decimal numbers
I wrote a shell bash script that can allow anyone to convert fractions into decimal numbers for downstream analyses. In my case, i was mainly interested on Major Allele (MAA) and Minor Allele (MIA) frequencies from the PoolSeq GWAS generated using Popoolation2.

How to run the command?
- Invoke `convert_fractions.sh` script for the job.

  
`./convert_fractions.sh -i path/to/the/input/file/input.txt -o path/to/the/input/file/output.txt -c ""`

Where
- -i defined the input file
- -o defined the output file
- -c defined the columns to process separated by space (eg. "5 8 10" for conversion of column 5, column 8 and column 10).

<!--
**Gadji-M/Gadji-M** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->
