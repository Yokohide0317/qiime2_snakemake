"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inputdir", "-i")

args = parser.parse_args()
inpDir = args.inputdir
"""

ONE = "01_data-import"
TWO = "02_adapter"

rule all:
    input: f"{ONE}/demux-paired-end.qzv"

rule manifest:
    input: "raw-data"
    output: "manifest.tsv"
    shell: "seqfu metadata --format manifest --split _S raw-data/ > manifest.tsv"


# {{{ =========== 01_IMPORT =============
rule import_a:
    input: "manifest.tsv"
    output: f"{ONE}/demux-paired-end.qza"
    shell: f"""
            qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {{input}} \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path {{output}}
            """

rule import_v:
    input: f"{ONE}/demux-paired-end.qza"
    output: f"{ONE}/demux-paired-end.qzv"
    shell: f"""
            qiime demux summarize \
            --i-data {{input}} \
            --o-visualization {{output}}
            """

# }}}

# {{{ =========== 02_ADAPTER =============
rule adapter_a:
    input: f"{ONE}/demux-paired-end.qza"
    output: f"{TWO}/trimmed-seqs.qza"
    shell: f"""
            qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {{input}} \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path {{output}}
            """

rule adapter_v:
    input: f"{ONE}/demux-paired-end.qza"
    output: f"{ONE}/demux-paired-end.qzv"
    shell: f"""
            qiime demux summarize \
            --i-data {{input}} \
            --o-visualization {{output}}
            """
# }}}