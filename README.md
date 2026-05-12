# `intom`: Perturbational robustness estimation using `int`egrative `om`ics data <a href="#"><img src="https://openmoji.org/data/color/svg/E25B.svg" align="right" width="20%"></a>

<!-- ![GitHub Repo stars](https://img.shields.io/github/stars/rraadd88/intom?style=for-the-badge)
<a href="">[![Downloads](https://img.shields.io/pypi/dm/intom?style=for-the-badge)](https://pepy.tech/project/intom)</a> -->
<!-- <a href="">[![DOI](https://img.shields.io/badge/DOI-link-blue?style=for-the-badge)](https://doi.org/10.5281/zenodo.20129325)</a> -->
<a href="">[![build](https://img.shields.io/github/actions/workflow/status/rraadd88/intom/build.yml?style=for-the-badge)](https://github.com/rraadd88/intom/actions/workflows/build.yml)</a>
<a href="">[![GNU Licence](https://img.shields.io/github/license/rraadd88/intom.svg?style=for-the-badge)](https://github.com/rraadd88/intom/blob/main/LICENSE)</a>
<a href="">[![Issues](https://img.shields.io/github/issues/rraadd88/intom.svg?style=for-the-badge)](https://github.com/rraadd88/intom/issues)</a>

<img width="781" height="306" alt="image" src="https://github.com/user-attachments/assets/c102c51f-057c-4c24-8dd2-516296ef89c9" />

## Usage

### Using perturbational data directly

    intom diff --input-path "inputs/sids.tsv" --output-path "outputs/diff_cli.yaml" --values-path "inputs/freq.tsv" --col-id "t id" --col-group "g id" --cols-value '["freq", "pert"]' --pms-diff '{"func": "ttest"}'
        
    ⚙️ Processing ..
    100%|████████████████|
    🗂️ Saving outputs ..
    ✅ Done.

<details>
    <summary><kbd>intom diff --help</kbd></summary>
    
    usage: intom diff [-h] [--input-path INPUT_PATH] [--output-path OUTPUT_PATH] [--col-id COL_ID] [--col-group COL_GROUP] [--col-sid COL_SID] [--values-path VALUES_PATH] [--cols-value COLS_VALUE] [-p PMS_DIFF] [-b BU_SCORE_GT] [--kws-scoring KWS_SCORING] [--kws-ro KWS_RO]
                      [--output-dir-path OUTPUT_DIR_PATH] [-t TABLE_EXT] [--verbose VERBOSE] [-r REF_PATH] [--ids-path IDS_PATH] [--mod-cfgs MOD_CFGS] [--mod-names MOD_NAMES] [-s SIDS] [--kws-scoring-flt KWS_SCORING_FLT] [--col-ids-skip-in COL_IDS_SKIP_IN]
                      [--kws-run-tasks-default KWS_RUN_TASKS_DEFAULT]
    
    Differential analysis to robustness probabilities.
    
    Notes:
        See the `intom/examples/inputs` directory for example inputs.
        Use json-string format for the non-string parameters.
        Parameters could be provided using a dictionary file: `--input-path params.yaml`
    
    options:
      -h, --help            show this help message and exit
      --input-path INPUT_PATH
                            -
      --output-path OUTPUT_PATH
                            -
      --col-id COL_ID       -
      --col-group COL_GROUP
                            -
      --col-sid COL_SID     'sample id'
      --values-path VALUES_PATH
                            -
      --cols-value COLS_VALUE
                            -
      -p PMS_DIFF, --pms-diff PMS_DIFF
                            {'col_diffid': None, 'func': 'limma', 'kws_diff': {'robust': True, 'trend': False, 'env_name': None}, 'col_score': 'LFC'}
      -b BU_SCORE_GT, --bu-score-gt BU_SCORE_GT
                            0
      --kws-scoring KWS_SCORING
                            {'value_min': 0.1}
      --kws-ro KWS_RO       {}
      --output-dir-path OUTPUT_DIR_PATH
                            -
      -t TABLE_EXT, --table-ext TABLE_EXT
                            'pqt'
      --verbose VERBOSE     -
      -r REF_PATH, --ref-path REF_PATH
                            -
      --ids-path IDS_PATH   -
      --mod-cfgs MOD_CFGS   -
      --mod-names MOD_NAMES
                            -
      -s SIDS, --sids SIDS  -
      --kws-scoring-flt KWS_SCORING_FLT
                            {}
      --col-ids-skip-in COL_IDS_SKIP_IN
                            'ids'
      --kws-run-tasks-default KWS_RUN_TASKS_DEFAULT
                            {'kernel': 'intom', 'force': False, 'test': False, 'cpus': 1}
    
</details>

### (Optional) Pre-processing

    intom pair --input-path "inputs/tids.tsv" --output-path "outputs/pair_cli.yaml" --col-id "t id" --col-group "g id" --mod-names '["freq", "pert"]' --values-path "inputs/mods.h5mu" --pairs-path "inputs/tids_sim.tsv" --drop-ids-path "inputs/drop_ids.tsv" --drop-ids-ref-path "inputs/drop_ids_ref.tsv" --kws-pre '{"corr_min": 0.1}'

<details>
    <summary><kbd>intom pair --help</kbd></summary>
    
    usage: intom pair [-h] [--input-path INPUT_PATH] [-o OUTPUT_PATH] [--pairs-path PAIRS_PATH] [-s SIDS] [--drop-ids-path DROP_IDS_PATH] [--drop-ids-ref-path DROP_IDS_REF_PATH] [-v VALUES_PATH] [--mod-names MOD_NAMES] [--mod-cfgs MOD_CFGS] [--idt IDT]
                      [--col-sid COL_SID] [--col-group COL_GROUP] [--col-id COL_ID] [--pms-diff PMS_DIFF] [-b BU_SCORE_GT] [--kws-pre KWS_PRE] [--kws-pair KWS_PAIR] [--kws-scoring KWS_SCORING] [--kws-scoring-flt KWS_SCORING_FLT] [--kws-ro KWS_RO]
                      [--table-ext TABLE_EXT] [-w WD_PATH] [--kernel KERNEL] [--cpus CPUS] [-f] [--test]
    
    (Optional) Pre-processing to pair IDs.
    
    Notes:
        See the `intom/examples/inputs` directory for example inputs.
        Use json-string format for the non-string parameters.
        Parameters could be provided using a dictionary file: `--input-path params.yaml`
    
    options:
      -h, --help            show this help message and exit
      --input-path INPUT_PATH
                            -
      -o OUTPUT_PATH, --output-path OUTPUT_PATH
                            -
      --pairs-path PAIRS_PATH
                            -
      -s SIDS, --sids SIDS  -
      --drop-ids-path DROP_IDS_PATH
                            -
      --drop-ids-ref-path DROP_IDS_REF_PATH
                            -
      -v VALUES_PATH, --values-path VALUES_PATH
                            -
      --mod-names MOD_NAMES
                            -
      --mod-cfgs MOD_CFGS   -
      --idt IDT             -
      --col-sid COL_SID     'sample id'
      --col-group COL_GROUP
                            -
      --col-id COL_ID       -
      --pms-diff PMS_DIFF   {'col_diffid': None, 'func': 'limma', 'kws_diff': {'robust': True, 'trend': False, 'env_name': None}}
      -b BU_SCORE_GT, --bu-score-gt BU_SCORE_GT
                            0
      --kws-pre KWS_PRE     {}
      --kws-pair KWS_PAIR   {}
      --kws-scoring KWS_SCORING
                            {'value_min': 0.1}
      --kws-scoring-flt KWS_SCORING_FLT
                            {}
      --kws-ro KWS_RO       {}
      --table-ext TABLE_EXT
                            'pqt'
      -w WD_PATH, --wd-path WD_PATH
                            -
      --kernel KERNEL       'intom'
      --cpus CPUS           1
      -f, --force           False
      --test                False    
</details>

## Installation

To use with limma (recommended), refer to [setup.sh](./setup.sh)

    cd intom
    uv sync
    uv pip install .

## Citation
    
    @article{dandage_2026a,
    title = {Recursive mutational robustness in cancer through intra- and intergenic compensation},
    }

<!--
    author       = {Dandage, Rohan  AND ..},    
    month        = ,
    year         = 2026,
    publisher    = {},
    doi          = {},
    url          = {https://doi.org/},
    -->
