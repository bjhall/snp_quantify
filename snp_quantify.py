#!/usr/bin/env python
import io
import os
import sys
import gzip
import subprocess
from collections import defaultdict, Counter
from pprint import pprint as pp

def main():
    input_vcf = sys.argv[1]
    input_bam = sys.argv[2]
    snps = create_bed_file_with_informative_snps(input_vcf)
    correction_factors = get_correction_factors(snps)
    quantities = quantify_bam(input_bam, correction_factors, snps)

    for sample_id, quantity in quantities.items():
        print(f"{sample_id}\t{quantity}")



def quantify_bam(bam, correction_factors, snps):
    proc = subprocess.Popen(["samtools", "mpileup", bam, "-l", "tmp.bed", "--output-BP", "--output-MQ"], stdout=subprocess.PIPE)
    counts = defaultdict(int)

    # The mess below is a primitive parser for mpileup output :)
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        chrom, pos, ref, dp, base_str, qual_str, basepos_str, mapq_str = line.strip().split("\t")
        if int(dp) >= 4:
            continue
        snp_data = snps[chrom][int(pos)]
        idx = 0
        num = ""
        skip_chars = 0
        in_number = False
        for base in base_str:
            if skip_chars > 0:
                skip_chars -= 1
                continue
            base = base.upper()
            if in_number:
                if base.isnumeric():
                    num += str(base)
                else:
                    skip_chars = int(num)
                    in_number = False
                continue
            if base == "^":
                skip_chars = 1
                continue
            if base in ["$", "!"]:
                continue
            if base in ["+", "-"]:
                in_number = True
                num = ""
                continue
            qual = ord(qual_str[idx]) - 33
            mapq = ord(mapq_str[idx]) - 33
            idx += 1

            # If base matches the SNP base and quality if good enough, count it
            if base == snp_data["snp_base"] and qual >= 36 and mapq >= 60:
                counts[snp_data["sample"]] += int(snp_data["weight"])

    # Correct for uneven number of informative SNPs in the samples
    return {sample_id: count * correction_factors[sample_id] for sample_id, count in counts.items()}



def get_correction_factors(snps):
    # Count the number of informative SNPs per sample (A)
    snp_counts = Counter(snp['sample'] for chrom in snps.values() for snp in chrom.values())

    # Number of SNPs per sample, if uniformly distributed (B = total SNPs / number of samples)
    uniform_n_snps = sum(snp_counts.values()) / len(snp_counts)

    # For each sample, calculate a correction factor according to B/A
    return {sample_id: uniform_n_snps / sample_n_snps for sample_id, sample_n_snps in snp_counts.items()}



def read_cached_bed(bed_path):
    print("Reading cached snps...", file=sys.stderr)
    snps = defaultdict(dict)
    with open(bed_path) as fh:
        for line in fh:
            chrom, start, end, name = line.strip().split("\t")
            sample, snp_base, weight = name.split(";")
            snps[chrom][int(end)] = {'sample':sample, 'snp_base':snp_base, 'weight':weight}
    return snps


def create_bed_file_with_informative_snps(vcf_path):
    out_bed_path = vcf_path + ".informative_snps.bed"

    # If cached SNP file exists, use it instead of recalculating everything
    if os.path.exists(out_bed_path):
        return read_cached_bed(out_bed_path)

    snps = defaultdict(dict)
    out = open(out_bed_path, "w")
    with gzip.open(sys.argv[1], "rb") as vcf:
        for line in vcf:

            # Skip comments
            if line.startswith(b"##"):
                continue

            # Store header columns
            if line.startswith(b"#"):
                header = line.decode("utf-8").strip().split("\t")
                continue

            parts = line.decode("utf-8").strip().split("\t")

            # Skip sex chromosomes
            if parts[0] in ["chrY", "chrX"]:
                continue

            # Skip indels
            if len(parts[3]) != 1 or len(parts[4]) != 1:
                continue

            # Count the number of 0/0s and 0/1&1/1s
            zero, non_zero = 0, 0
            col_with_snp, zygosity = None, None
            for i in range(9, len(parts)):
                genotype = parts[i].split(":")[0].replace("|", "/")
                if genotype in ["0/1", "1/1"]:
                    non_zero += 1
                    col_with_snp = i
                    weight = 2 if genotype == "0/1" else 1 # Set weight of SNP to 2 if heterozygotic
                elif genotype == "0/0":
                    zero += 1

            # If exacly one sampled had the variant, and all others were 0/0, keep it!
            if non_zero == 1 and zero == len(parts)-10:

                # Write it to the cache file
                out.write("\t".join([
                    parts[0],
                    str(int(parts[1])-1),
                    parts[1],
                    header[col_with_snp] + ";" + parts[4] + ";" + str(weight),
                ]) + "\n")

                # Store it in a dict
                snps[parts[0]][int(parts[1])] = {'sample':header[col_with_snp], 'snp_base': parts[4], 'weight': weight}

    out.close()
    return snps

if __name__ == "__main__":
    main()
