#!/usr/bin/env python
import io
import os
import sys
import pandas as pd
import random
import subprocess

vcf = "cohort.gatk_vcf.filt500.sort.chr.vcf.gz"

samples = {
    'WH-3759-H9.md.subsampled.bam': 846879,
    'WH-3759-HS1001.md.subsampled.bam': 906706,
    'WH-3759-RC17.md.subsampled.bam': 938444,
}


dfs = []
for tot_reads in [300_000, 100_000, 10_000, 1_000]:
    for rep in range(0,50):
        fractions = list(random.sample(range(0, 100000), 3))
        fractions = [e/sum(fractions) for e in fractions]
        print(f"Processing: {tot_reads} reads, rep {rep+1} {fractions[0]}", file=sys.stderr)

        first = True
        if os.path.exists("tmp.sam"):
            os.remove("tmp.sam")

        # Merge reads from the three bams
        merged_out_sam = open("tmp.sam", "a")
        metadata = {}
        for bam, total_bam_reads in samples.items():
            sample_name = bam.split(".")[0]
            final_sample_fraction = fractions.pop()
            metadata[sample_name] = {'expected_fraction': final_sample_fraction}

            bam_fraction = (tot_reads * final_sample_fraction) / total_bam_reads
            cmd = ["samtools", "view", "-s", str(bam_fraction), bam]
            if first:
                cmd.append("-h")
            subprocess.run(cmd, stdout=merged_out_sam, stderr=subprocess.DEVNULL)
            first = False

        merged_out_sam.close()
        subprocess.run(["samtools", "view", "-b", "tmp.sam", "-o", "tmp.bam"])
        subprocess.run(["samtools", "sort", "tmp.bam", "-o", "tmp.sort.bam"])
        os.remove("tmp.sam")
        os.remove("tmp.bam")

        # Quantify with snp_quantify
        out = subprocess.check_output(["./snp_quantify.py", vcf,  "tmp.sort.bam"], stderr=subprocess.DEVNULL).decode("utf-8")
        df = pd.read_csv(io.StringIO(out), sep="\t")
        df['expected_fraction'] = df['sample'].map(lambda x: metadata[x]['expected_fraction'])
        df['observed_fraction'] = df['count'] / df['count'].sum()
        df['total_reads'] = tot_reads
        dfs.append(df)
        os.remove("tmp.sort.bam")

df = pd.concat(dfs)
df.to_csv("simulated.tsv", sep="\t")
