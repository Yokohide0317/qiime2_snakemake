# qiime2_snakemake

## Prepare
```:bash
conda activate qiime2-2022.11
mamba install seqfu
mamba install snakemake

# Download taxonomy database
wget https://data.qiime2.org/2023.2/common/silva-138-99-nb-classifier.qza
```

## config.yaml
| Key | Desc | Value (Default) |
| :---: | :---: | :---: |
| FWD | Forwardのアダプター配列 | String |
| REV | Reverseのアダプター配列 | String |
| CUT_COUNT_PERCENT | リードの後半部分で、Maxリード数よりもリード数が少なくなった時に切る | Float (0.9) |
| CUT_QUALITY_PERCENT | QCで目安にするクォリティ値 | Int (25) |
| CUT_QUALITY_LOCATION | QCで目安にする部位 | Choice(25, 50, 75) (50) |
| TAXONOMY_DATABASE | 学習済みQiime2用分類器 | Path ('./silva-138-99-nb-classifier.qza') |
| CUTADAPT_CORES | cutadaptに使用するcore数 | Int (4) |
| DADA2_CORES | DADA2に使用するcore数 | Int (4) |
| TAXONOMY_CORES | taxonomyに使用するcore数 | Int (4) |

## Exec
```:bash
mkdir raw-data
mv <SAMPLEPATH/*_fastq.gz> raw-data/
snakemake --cores N
```

## Adapter Examples
```
# V1-V2
FWD: "AGRGTTTGATYMTGGCTCAG"
REV: "TGCTGCCTCCCGTAGGAGT"

# V3-V4
FWD: "CCTACGGGNGGCWGCAG"
REV: "GACTACHVGGGTATCTAATCC"
```