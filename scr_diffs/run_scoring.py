#!/usr/bin/env python
# coding: utf-8
# # Scoring with ref/similar samples are known
#
# ```mermaid
# graph LR;
#     bvb["backup vs backup"]
#         --> bvnl["backup vs no-loss"] & bvl["backup vs loss"]
#     bvnl
#         --> nl_in_nm["expression of\nthe transcripts\nin no loss"] & b_in_nm["expression of\nthe transcripts\nin backup"]
#     bvl
#         --> nl_in_m["expression of\nthe transcripts\nin no loss"] & b_in_m["expression of\nthe transcripts\nin backup"]
# ```
#
# Notes: refactored for reusing where raw inputs are not available.
#
# Input formats:
#
#
# **1 file:**
#
# - `input_path`: Merged data
#     e.g.
#
#         comparison	gene id	ids	ids var2	mRNA expression (TPM)	sample id	transcript id
#
# **2 files:**
#
# - `input_path`: Test data
#     e.g.
#
#         gene id	ids	ids var2	mRNA expression (TPM)	sample id	transcript id
#
# - `ref_path`: with ref data with data for similar samples provided
#
#        gene id	ids	mRNA expression (TPM)	rp	sample id	sample id deleteriousness	transcript id
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    col_id=None,
    col_group=None,
    col_sid="sample id",
    values_path=None,
    cols_value=None,
    pms_diff={
        "col_diffid": None,
        "func": "limma",
        "kws_diff": {"robust": True, "trend": False, "env_name": None},
        "col_score": "LFC",
    },
    bu_score_gt=0,
    kws_scoring={"value_min": 0.1},
    kws_ro={},
    output_dir_path=None,
    table_ext="pqt",
    verbose=None,
    ref_path=None,
    ids_path=None,
    mod_cfgs=None,
    mod_names=None,
    sids=None,
    kws_scoring_flt={},
    col_ids_test=None,
    col_ids_skip_in="ids",
    kws_run_tasks_default={"kernel": "intom", "force": False, "test": False, "cpus": 1},
):
    if input_path is not None and output_path is None:
        ## should be params
        from roux.workflow.task import pre_params

        params = pre_params(input_path)
        # if len(params)!=1:
        #     logging.info(f"processing {len(params)} pms")
        from roux.workflow.function import call_with_kws

        for i, _kws in enumerate(params):
            call_with_kws(
                run,
                _kws,
                **{
                    k: v
                    for k, v in dict(
                        verbose=locals().get("verbose"),
                        force=locals().get("force"),
                    ).items()
                    if v is not None
                },
            )
        return

    ## for cli
    assert input_path is not None, "check docs using --help"
    assert output_path is not None

    import logging

    logging.basicConfig(
        level=verbose
        if isinstance(verbose, str)
        else "INFO"
        if verbose == True
        else 29,
        force=True,
    )
    from pathlib import Path
    from roux.lib.io import to_dict, read_table, to_table

    ## workflow functions from roux
    from roux.workflow.task import run_tasks  # noqa
    from roux.workflow.log import test_params
    from roux.workflow.io import check_for_exit

    ## for json formatted pms from CLI
    def set_type(pm):
        import json

        if isinstance(pm, str):
            pm = json.loads(pm)
            assert not isinstance(pm, str), pm
        return pm

    cols_value = set_type(cols_value)
    pms_diff = set_type(pms_diff)
    kws_scoring = set_type(kws_scoring)
    kws_run_tasks_default = set_type(kws_run_tasks_default)

    ## source dir.
    from intom.run import get_src_path

    _src_dir_path = get_src_path()

    if pms_diff["func"] == "limma":
        pms_diff["kws_diff"]["env_name"] = pms_diff["kws_diff"].get(
            "env_name", kws_run_tasks_default["kernel"]
        )
        from intom.stats.diffs import get_ro

        get_ro(
            pms_diff["kws_diff"]["env_name"],
            verbose=True,
        )

    if output_dir_path is None:
        output_dir_path = Path(output_path).with_suffix("").as_posix()
    logging.info(f"Out. dir. path: {output_dir_path}")

    if col_ids_test is None:
        col_ids_test = f"{col_ids_skip_in} {cols_value[1]}"

    ## diff. by
    col_diffid = pms_diff.get("col_diffid", col_group)
    pms_diff["col_diffid"] = col_diffid

    cols_id = dict(
        col_group=col_group,
        col_id=col_id,
        col_sid=col_sid,
        col_ids_test=col_ids_test,
    )
    assert all([(v is not None) for k, v in cols_id.items()]), [
        (k, (v is not None)) for k, v in cols_id.items()
    ]
    # ## Inputs

    cfgt_path = """
    ## backup upreg.    
    bu:
        dir_path: ${..dir_path}/01_backup_upregulation
        path: ${.dir_path}.${table_ext}
    # compensation by backup
    bc:
        dir_path: ${..dir_path}/02_backup_compensation
        path: ${.dir_path}.${table_ext}
    ## robustness, combined output
    ro:
        dir_path: ${..dir_path}/03_robust
        path: ${.dir_path}.${table_ext}
    """
    from roux.workflow.io import read_config  # noqa

    cfg = read_config(
        cfgt_path,
        inputs=dict(
            dir_path=output_dir_path,
            table_ext=table_ext,
        ),
    )
    cfg
    # ### Comparisons

    df0_sids = read_table(input_path)

    mappedby_values = cols_value[0] in df0_sids
    if mappedby_values:
        comps_mappedby_values_path = input_path
        del df0_sids
    else:
        ## need to map the data
        comps_mappedby_values_path = (
            f"{output_dir_path}/01_backup_upregulation/00_combined_diff_input_raw.pqt"
        )
        mappedby_values = Path(comps_mappedby_values_path).exists()
    # #### Mapping values to the comparisons

    if not mappedby_values:
        ## explode ids
        if col_ids_skip_in in df0_sids and col_id not in df0_sids:
            df1_sids = (
                df0_sids.assign(
                    **{
                        col_id: lambda df: df[col_ids_skip_in].str.split(";"),
                    }
                )
                .explode(col_id)
                .log.head()
            )
        else:
            df1_sids = df0_sids.copy()

    if not mappedby_values:
        ids = df1_sids[col_id].unique().tolist()
        sids = df1_sids[col_sid].unique().tolist()
        assert len(ids) > 0, len(ids)
        assert len(sids) > 0, len(sids)

    if not mappedby_values:
        from intom.omics.io_prep import read_table_prep

        # %run ../modules/omics/io_prep.py
        df0_values = read_table_prep(
            values_path,
            col_group=col_group,
            col_id=col_id,  # col_id_ext if col_id_ext is not None else col_id,
            mod_name=cols_value[0],
            sids=sids,
            ids=ids,
        )
    # #### Mapping values
    #
    # Format:
    #
    #     comparison    gene id    ids    ids var2    sample id    transcript id   mRNA expression (TPM)

    if not mappedby_values:
        df1_comps = df1_sids.log.merge(
            right=df0_values.drop(["comparison", col_group], axis=1, errors="ignore"),
            how="inner",
            on=[col_id, col_sid],
            validate="m:1",
        ).log.head()

    if not mappedby_values:
        df1_comps = df1_comps.rd.assert_dense(
            # .check_dups(
            subset=df1_comps.drop(
                [cols_value[0]], axis=1, errors="ignore"
            ).columns.tolist()
        )

    if not mappedby_values:
        comps_mappedby_values_path = to_table(
            df1_comps,
            comps_mappedby_values_path,
        )
        del df0_values
    # ## Backup-upregulation

    params = []
    params.append(
        {
            **dict(
                input_path=comps_mappedby_values_path,  # cfg["inputs"]['02_filteredby_pairs']['path'],
                output_path=cfg["bu"]["path"],
                ref_path=ref_path,  # cfg["inputs"]['03_backup_vs_noloss']['path'],
                cols=cols_value,
                values_path=values_path,
                col_diffid=col_diffid,
            ),
            **kws_scoring,
            **pms_diff,
            **cols_id,
        }
    )
    test_params(params)

    outputs = run_tasks(
        f"{_src_dir_path}/scr_diffs//01_upregulation.py",
        params=params,
        **kws_run_tasks_default,
    )
    outputs

    check_for_exit(
        params[0]["output_path"],
        outp=output_path,
        data=cfg,
    )
    # ### Validation: filtering to significant effects and validations by trace-back

    params = []

    params.append(
        {
            **dict(
                input_path=cfg["bu"]["path"],
                output_path=f"{Path(cfg['bu']['path']).with_suffix('').as_posix()}/02_diff_top.{table_ext}",
                ## for trace-back
                values_path=f"{Path(cfg['bu']['path']).parent.as_posix()}/00_combined.pqt",
                cols=cols_value,
                src_dir_path=_src_dir_path,
            ),
            **cols_id,
            **kws_scoring,
            **kws_scoring_flt,
            **{k: pms_diff[k] for k in ["col_score"] if k in pms_diff},
        }
    )
    test_params(params)

    try:
        outputs = run_tasks(
            f"{_src_dir_path}/scr_diffs//011_filtering_to_signi_effects.py",
            params=params,
            **kws_run_tasks_default,
        )
    except Exception as e:
        logging.warning(e)
    # ## Compensation (acts as an independent validation too)

    params = []
    params.append(
        {
            **dict(
                input_path=cfg["bu"]["path"],
                output_path=cfg["bc"]["path"],
                exp_dir_path=(
                    mod_cfgs[mod_names[0]]["path"]
                    if values_path is None
                    else values_path
                ),
                ids_path=ids_path,
                sids=sids,
                cols_value=cols_value,
                score_gt=bu_score_gt,
            ),
            **cols_id,
            **{k: pms_diff[k] for k in ["col_score"] if k in pms_diff},
        }
    )
    test_params(params)

    # %%time
    run_tasks(
        f"{_src_dir_path}/scr_diffs//02_compensation.py",
        params=params,
        **kws_run_tasks_default,
    )

    check_for_exit(
        params[0]["output_path"],
        outp=output_path,
        data=cfg,
    )
    # ## Robustnness

    params = []
    params.append(
        {
            **dict(
                input_path=cfg["bc"]["path"],
                output_path=cfg["ro"]["path"],
                cols_value=cols_value,
            ),
            **kws_ro,
            **cols_id,
            **{k: pms_diff[k] for k in ["col_score"] if k in pms_diff},
        },
    )
    test_params(params)

    # %%time
    run_tasks(
        f"{_src_dir_path}/scr_diffs//03_robust.py",
        params=params,
        **kws_run_tasks_default,
    )

    check_for_exit(
        params[0]["output_path"],
        outp=output_path,
        data=cfg,
    )
    # ## Outputs

    to_dict(
        cfg,
        output_path,
    )
    return output_path


## CLI-setup
parser = argh.ArghParser()
parser.add_commands(
    [
        run,
    ]
)
if __name__ == "__main__":  # and sys.stdin.isatty():
    parser.dispatch()
