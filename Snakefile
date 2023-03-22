import pandas as pd


configfile: "config.yaml"

ONE = "01_data-import"
TWO = "02_adapter"
THR = "03_denoise"

translateDict = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "m": "k",
    "r": "y",
    "y": "r",
    "k": "m",
    "v": "b",
    "h": "d",
    "d": "h",
    "b": "v",
}
FWD: str = config["FWD"].strip()
REV: str = config["REV"].strip()
FWDRC: str = ''.join(list(reversed(FWD.translate(str.maketrans(translateDict)))))
REVRC: str = ''.join(list(reversed(REV.translate(str.maketrans(translateDict)))))

print(FWDRC)
print(REVRC)

CUT_COUNT_PERCENT = config["CUT_COUNT_PERCENT"]
CUT_QUALITY_PERCENT = config["CUT_QUALITY_PERCENT"]
CUT_QUALITY_LOCATION = config["CUT_QUALITY_LOCATION"]

# QCの数値出し
def cut_edge(_df: pd.DataFrame, _CUT_COUNT_PERCENT: float):
    max_count: int = int(_df.loc["count"][0])

    drop_colnum: int = int(_df.columns[-1])
    for col, data in enumerate(_df.loc["count"]):
        if data < max_count * _CUT_COUNT_PERCENT:
            drop_colnum = int(col)
            break

    filt_df = _df.drop([str(col) for col in _df.columns[drop_colnum:]], axis=1)
    return filt_df

def quality_check(_df: pd.DataFrame, _CUT_QUALITY_PERCENT: int, _CUT_QUALITY_LOCATION: int) -> int:
    target_row = _df.loc[f"{str(_CUT_QUALITY_LOCATION)}%"]

    trunc_loc: int = int(target_row.keys()[-1])
    for col, data in target_row.items():
        if data < _CUT_QUALITY_PERCENT:
            trunc_loc = int(col) -1
            break

    return trunc_loc



rule all:
    input: f"{THR}/stats-dada2.qza"

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
            qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {{input}} \
            --p-front-f {FWD} \
            --p-front-r {REV} \
            --p-adapter-f {FWDRC} \
            --p-adapter-r {REVRC} \
            --p-times 10 \
            --p-match-read-wildcards \
            --p-match-adapter-wildcards \
            --p-minimum-length 100 \
            --p-discard-untrimmed \
            --o-trimmed-sequences {{output}}
            """

rule adapter_v:
    input: f"{TWO}/trimmed-seqs.qza"
    output: f"{TWO}/trimmed-seqs.qzv"
    shell: f"""
            qiime demux summarize \
            --i-data {{input}} \
            --o-visualization {{output}}
            """

rule adapter_o:
    input: f"{TWO}/trimmed-seqs.qzv"
    output: f"{TWO}/output/forward-seven-number-summaries.tsv"
    shell: f"""
            qiime tools export \
            --input-path {{input}} \
            --output-path {TWO}/output
            """
# }}}

# {{{ =========== 03_denoise =============
rule denoise_a:
    input: f"{TWO}/trimmed-seqs.qza"
    output: f"{THR}/table-dada2.qza", f"{THR}/rep-seqs-dada2.qza", f"{THR}/stats-dada2.qza"
    params:
        #tsv_path_f = f"{TWO}/output/forward-seven-number-summaries.tsv"
        trunc_len = [
            quality_check(cut_edge(pd.read_csv(f"{TWO}/output/forward-seven-number-summaries.tsv", sep="\t", index_col=0), _CUT_COUNT_PERCENT=CUT_COUNT_PERCENT), _CUT_QUALITY_PERCENT=CUT_QUALITY_PERCENT, _CUT_QUALITY_LOCATION=CUT_QUALITY_LOCATION),
            quality_check(cut_edge(pd.read_csv(f"{TWO}/output/reverse-seven-number-summaries.tsv", sep="\t", index_col=0), _CUT_COUNT_PERCENT=CUT_COUNT_PERCENT), _CUT_QUALITY_PERCENT=CUT_QUALITY_PERCENT, _CUT_QUALITY_LOCATION=CUT_QUALITY_LOCATION)
        ]
        #tsv_path_r = f"{TWO}/output/reverse-seven-number-summaries.tsv"
        #trunc_len_r = quality_check(cut_edge(pd.read_csv(tsv_path_r, sep="\t", index_col=0), CUT_COUNT_PERCENT), CUT_QUALITY_PERCENT, CUT_QUALITY_LOCATION)

    shell: f"""
            qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {{input}} \
            --p-trunc-len-f {{params.trunc_len[0]}} \
            --p-trunc-len-r {{params.trunc_len[1]}} \
            --o-table {THR}/table-dada2.qza \
            --o-representative-sequences {THR}/rep-seqs-dada2.qza \
            --o-denoising-stats {THR}/stats-dada2.qza
            """
# }}}