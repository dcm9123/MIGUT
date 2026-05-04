###April 19th, 2026

After waiting for months, I finally have the Parkinson's data with me. I will be presenting some basic results to the public this upcoming Friday at the PAA (Lunch & Learn). The directory and samples I am working with are located in ARC, specifically in `/bulk/IMCshared/daniel/novaseq/data`, where I have a total of 186 samples divided into two .fastq files (372 .fastq.gz files in total).

I tried working with the existing metaphlan4 pipeline from the Sycuro Lab, but I am having compatibility issues with my Snakemake version vs theirs. Because of that, I decided to run my own analysis in my own semi-automatic workflow by using more updated software and database versions. The first thing I noticed is that the `config` file for Snakemake had a very short adapter sequence (short for Illumina, at least). After visualizing the .fastqc and .multiqc reports of the data, I could see that whomever run that pipeline, did wrong the use adapter sequences, as multiqc was showing a very-high level of adapter and representative sequences that should not be there after running cutadapt and prinseq. Because of that, I decided to look for the right ones and start over with my own results.

###April 27th, 2026

I had a chance to present my research at PAA, the talk went well. I had to rush some analyses and data generation to do so. Here is what I did in a nutshell:

1. With Kevin's help, I was able to fix the Snakemake `metqc` pipeline and run it on the data. I had to change the adapter sequences MULTIPLE TIMES, as we couldn't get the right ones. We kneew we were having adapter contamination as both the `multiqc` and individual `fastqc` reports were showing a very high level of adapter sequences. After trying multiple adapter sequences, we finally found the right ones and got the contamination to acceptable levels. We also had to add the `cutadapt` option to trim multiple GGGGGG homopolymers that seem to be an error with NextSeq and NovaSeq instruments. The folder with the latest results (Apr 23rd, 2026) is located in `/bulk/IMCshared_bulk/daniel/parkinsons/novaseq/metqc_output_daniel_snakemake`. This pipeline ran `cutadapt`, `prinseq`, `fastqc`, `multiqc`, and `bmtagger` to (i) remove adapter sequences and polyG homopolymers, (ii) filter out low-quality reads, and (iii) remove human contamination. We looked at the fastqc results, and the quality seemed much better.

2. The next step was to run `metaphlan4` on the cleaned, high-quality, non-human reads to get taxonomic profiles. The folder with the latest results. For that, I used a script located in `/bulk/IMCshared_bulk/daniel/parkinsons/novaseq/metqc_output_daniel_snakemake/utils/one_script_persample_metaphlan4.sh` that runs `metaphlan4`. I am also moving the scripts that automate this script submission into a folder in this repository. Its called `toolbox`. 


### May 4th, 2026
Whilst preparing my QC report on my reads, I noticed that the GC content of these files was relatively high, high enough to trigger a warning on the `multiqc` report. After noticing that, I decided to look into the forward and reverse .fastq files to see if there were any polyGs still present: `for file in *; do grep -c "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" ${file}; done;`. This gave me the following results:

```
10
13330
5
5940
11
4082
0
72
8
7397
3
10551
```

Every odd line represents the forward .fastq file, and every even line represents the reverse .fastq file. As you can see, there are still a lot of polyGs in the data, which is likely the reason for the high GC content. I tried adding a new parameter to `cutadapt`, specifically `-e 0.0 -G{20}`, which should trim any homopolymer of Gs that is 20 or more bases long. However, this did not seem to work, as I still had a lot of polyGs in the data. I decided to switch to `bbduk`, where I will run a two-step approach by first removing any adapters, and second, removing any polyGs. This is the command I will try first:

`for file in *_1.fastq.gz; do bbduk.sh in1=${file} in2=${file%%_1.fastq.gz}_2.fastq.gz out1=../bbduk_testing/{file%%_1.fastq.gz}_1_bbduk_filt.fastq.gz out2=../bbduk_testing/${file%%_1.fastq.gz}_2_bbduk_filt.fastq.gz trimq=20 qtrim=rl minlen=60 maq=20 maxns=2 hdist=1 k=27 stats=../bbduk_testing/${file%%_1.fastq.gz}_stats.txt tbo tpe threads=32 ref=adapters; done; `


