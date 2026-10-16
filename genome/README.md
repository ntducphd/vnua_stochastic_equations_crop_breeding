# Genome measurements

The *a1*–*sh2* distances (Section 3.3) and the chromosome 4 diversity contrast around *su1* (Section 4.4)
are read from the tables in `results/tables/`:

    python -I -B genome/scripts/genome_numbers.py

The script prints each quoted value with its definition and writes them to `results/genome_numbers.txt`.

| File | Contents |
|---|---|
| `scripts/07_population_differentiation.py` | Hudson F<sub>ST</sub> and expected heterozygosity in 1-Mb windows from the GBS calls of Hu et al. (2021) |
| `scripts/13_a1_sh2_interval.py` | gene order, orientation and the *a1*–*sh2* distance in 17 assemblies, with the assembly gaps in the interval |
| `scripts/01_download_inputs.py` | downloads the inputs of 07 and 13 and checks them against `data/MANIFEST_inputs_md5.tsv` |
| `scripts/genome_numbers.py` | the values quoted in the paper |
| `data/curated/` | key genes, assemblies and founder lines, and the F<sub>2</sub> kernel counts of Ravichandran et al. (2026) Table S1 |
| `results/tables/` | window statistics (B73 RefGen_v4), key-gene positions (v4 and v5) and the *a1*–*sh2* interval per assembly |

## Sources

- Hu et al. (2021) *Nat. Commun.* 12:1227 (CC BY): GBS calls `noSeed_denovo_publicSweet.vcf.gz` (CyVerse
  dataset `UF_resende_Sweet_corn_2020`) and the accession groups of Supplementary Data 4. Matched to the
  calls, the groups give 798 sweet corn and 520 temperate field corn accessions and 859,631 SNPs.
- MaizeGDB (https://download.maizegdb.org, accessed 12 October 2026): gene maps, gene tables, annotations
  and genome sequences of the assemblies.

## Rebuilding the tables

1. `python -I genome/scripts/01_download_inputs.py --local <folder>` downloads about 11.7 GB into `data/raw/`.
2. `python -I genome/scripts/07_population_differentiation.py`
3. `python -I genome/scripts/13_a1_sh2_interval.py` also needs the *sh2* k-mer anchor tables
   `kmer_hits_near_locus.csv.gz` and `parents_kmer_hits_near_locus.csv.gz` in `results/tables/`, produced by
   scanning the genome sequences; they are not included.
