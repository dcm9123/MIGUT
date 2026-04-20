April 19th, 2026

After waiting for months, I finally have the Parkinson's data with me. I will be presenting some basic results to the public this upcoming Friday at the PAA (Lunch & Learn). The directory and samples I am working with are located in ARC, specifically in `/bulk/IMCshared/daniel/novaseq/data`, where I have a total of 186 samples divided into two .fastq files (372 .fastq.gz files in total).

I tried working with the existing metaphlan4 pipeline from the Sycuro Lab, but I am having compatibility issues with my Snakemake version vs theirs. Because of that, I decided to run my own analysis in my own semi-automatic workflow by using more updated software and database versions. The first thing I noticed is that the `config` file for Snakemake had a very short adapter sequence (short for Illumina, at least). After visualizing the .fastqc and .multiqc reports of the data, I could see that whomever run that pipeline, did wrong the use adapter sequences, as multiqc was showing a very-high level of adapter and representative sequences that should not be there after running cutadapt and prinseq. Because of that, I decided to look for the right ones and start over with my own results.

 
