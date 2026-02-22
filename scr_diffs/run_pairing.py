#!/usr/bin/env python
# coding: utf-8
# # Pairing.
#
# Note: This script processes one group of samples e.g. one cancer subtype.
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    pairs_path=None,
    sids=None,
    drop_ids_path=None,
    drop_ids_ref_path=None,
    values_path=None,
    mod_names=None,
    mod_cfgs=None,
    idt=None,
    col_sid="sample id",
    col_group=None,
    col_id=None,
    pms_diff={
        "col_diffid": None,
        "func": "limma",
        "kws_diff": {"robust": True, "trend": False, "env_name": None},
    },
    bu_score_gt=0,
    kws_pre={},
    kws_pair={},
    kws_scoring={"value_min": 0.1},
    kws_scoring_flt={},
    kws_ro={},
    table_ext="pqt",
    wd_path=None,
    kernel="intom",
    cpus=1,
    verbose=None,
    force=False,
    test=False,
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

    ## logging functions
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

    ## data functions
    ## system functions from roux
    from roux.lib.sys import read_ps
    from roux.lib.io import read_table, to_dict  # , to_table

    ## visualization functions from roux
    ## workflow functions from roux
    from roux.workflow.task import run_tasks  # noqa
    from roux.workflow.log import test_params
    from roux.workflow.io import check_for_exit
    import nest_asyncio

    nest_asyncio.apply()

    ## source dir.
    from intom.run import get_src_path

    _src_dir_path = get_src_path()
    # ## Inputs

    ## validating parameters
    assert "kws_diff" not in locals(), "provide diff pms in pms_diff"
    assert output_path.endswith((".yml", ".yaml")), output_path

    if wd_path is not None:
        import os

        os.chdir(wd_path)
        print(wd_path)

    # logging.info(f"Out. dir. path: {output_dir_path}")
    from roux.workflow.io import set_outputs

    output_dir_path, output_paths = set_outputs(
        output_path=locals().get("output_path"),
        output_paths=locals().get("output_paths"),
    )
    # ### Modalities

    if values_path is None:
        if isinstance(mod_cfgs, (str, list)):
            mod_cfgs = read_ps(mod_cfgs)
            mod_cfgs = {
                str(Path(p).stem): dict(path=p, col=str(Path(p).stem)) for p in mod_cfgs
            }
    else:
        import mudata as md

        dt1 = md.read(values_path, backed=True)
        logging.info(dt1)
        if values_path.endswith(".h5ad"):
            mod_cfgs = dt1.uns
        else:
            mod_cfgs = {
                k: {k_: v_ for k_, v_ in dt1.mod[k].uns.items() if k_ in ["col"]}
                for k in dt1.mod_names
            }
        dt1.file.close()
    logging.info(mod_cfgs)

    if isinstance(mod_names, str):
        mod_names = mod_names.split("--")
    assert len(mod_names) == 2, mod_names
    logging.info(mod_names)

    if mod_names[0] not in mod_cfgs:
        mod_names = [
            [k for k, d in mod_cfgs.items() if d.get("col") == c][0] for c in mod_names
        ]
        logging.info(f"adjusted: mod_names={mod_names}")
    # #### Columns

    if idt is None:
        ## family members
        df00 = read_table(input_path)
        print(df00.head(1))
        idt = col_id.split(" id")[0]
        print(idt)
        assert not df00[[col_id, col_group]].duplicated().any()

    pms_cols_id = dict(
        col_group=col_group,
        col_id=col_id,
        col_sid=col_sid,
    )

    cols_value = [mod_cfgs[s].get("col", s) for i, s in enumerate(mod_names)]
    logging.info(cols_value)

    kws_run_tasks_default = dict(
        kernel=kernel,
        force=force,
        test=test,
        cpus=cpus,
    )
    # ### Paths

    cfgt_path = """
    dir_path: #from input 
    # file table_extension
    table_ext: pqt
    inputs:
        dir_path: ${..dir_path}/00_inputs
        00_vars_filtered:
            dir_path: ${..dir_path}/00_vars_filtered
            path: ${.dir_path}.${table_ext}
        01_ids_paired:
            dir_path: ${..dir_path}/01_ids_paired
            path: ${.dir_path}.${table_ext}
        02_filteredby_pairs:
            dir_path: ${..dir_path}/02_filteredby_pairs
            path: ${.dir_path}.${table_ext}
        03_backup_vs_noloss:
            dir_path: ${..dir_path}/03_backup_vs_noloss
            path: ${.dir_path}.${table_ext}
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
    from roux.workflow.io import read_config

    cfg = read_config(
        cfgt_path,
        inputs=dict(
            dir_path=output_dir_path,
            table_ext=table_ext,
        ),
    )
    cfg

    if isinstance(sids, str):
        sids = read_table(sids)["sample id"].tolist()
    # ## Input for scoring
    # ### Scoring inputs with paired values

    params = []
    params.append(
        {
            **dict(
                input_path=input_path,
                output_path=cfg["inputs"]["00_vars_filtered"]["path"],
                mod_names=mod_names,
                ## todo: provide mudata instead
                input_paths=(
                    values_path
                    if values_path is not None
                    ## for back-compatibility
                    else {
                        **{
                            f"var{i+1}": mod_cfgs[var]["path"]
                            for i, var in enumerate(mod_names)
                        },
                    }
                ),
                sids=sids,
                cols=cols_value,
                #     None if values_path is not None else
                #     ## for back-compatibility
                #     cols_value
                # ),
                ## for diffs
                paired=False,
                col_group=col_group,
                cols_id=[col_id],
                ## filtering
                drop_ids_path=drop_ids_path,
            ),
            **kws_pair,
        }
    )
    test_params(params)

    run_tasks(
        f"{_src_dir_path}/scr_scripts//merge_variables.py",
        params=params,
        **kws_run_tasks_default,
    )

    check_for_exit(
        params[0]["output_path"],
        outp=output_path,
        data=cfg,
    )
    # ### Pairing backup to loss ids

    params = []
    params.append(
        {
            **dict(
                input_path=cfg["inputs"]["00_vars_filtered"]["path"],
                output_path=cfg["inputs"]["01_ids_paired"]["path"],
                cols=cols_value,
            ),
            **pms_cols_id,
        }
    )
    test_params(params)

    run_tasks(
        f"{_src_dir_path}/scr_diffs//000_pairing_backup_and_loss_ids.py",
        params=params,
        **kws_run_tasks_default,
    )

    check_for_exit(
        params[0]["output_path"],
        outp=output_path,
        data=cfg,
    )
    # ### Filtering id pairs by the search space of interacting ids

    params = []
    params.append(
        {
            **dict(
                input_path=cfg["inputs"]["01_ids_paired"]["path"],
                output_path=cfg["inputs"]["02_filteredby_pairs"]["path"],
                pairs_path=pairs_path,
                cols=cols_value,
            ),
            **pms_cols_id,
        }
    )
    test_params(params)

    run_tasks(
        f"{_src_dir_path}/scr_diffs//001_filtering_pairs.py",
        params=params,
        **kws_run_tasks_default,
    )

    check_for_exit(
        params[0]["output_path"],
        outp=output_path,
        data=cfg,
    )
    # ### Backup vs no-loss sample and values

    params = []
    params.append(
        {
            **dict(
                input_path=cfg["inputs"]["02_filteredby_pairs"]["path"],
                output_path=cfg["inputs"]["03_backup_vs_noloss"]["path"],
                cols=cols_value,
                values_path=(
                    values_path
                    if values_path is not None
                    ## for back-compatibility
                    else mod_cfgs[mod_names[0]]["path"]
                ),
                sids=sids,
                values_flt_path=cfg["inputs"]["00_vars_filtered"]["path"],
                ids_path=input_path,
                drop_ids_path=drop_ids_ref_path,
            ),
            **pms_cols_id,
            **kws_pre,
        }
    )
    test_params(params)

    run_tasks(
        f"{_src_dir_path}/scr_diffs//002_backup_vs_no_loss_samples.py",
        params=params,
        **kws_run_tasks_default,
    )

    check_for_exit(
        params[0]["output_path"],
        outp=output_path,
        data=cfg,
    )
    # ## Scoring

    outp = Path(cfg["ro"]["path"]).with_suffix(".yaml").as_posix()
    params = []
    params.append(
        dict(
            ## params
            input_path=cfg["inputs"]["02_filteredby_pairs"]["path"],
            output_path=outp,
            ids_path=input_path,
            col_id=pms_cols_id["col_id"],
            col_group=pms_cols_id["col_group"],
            col_sid=pms_cols_id["col_sid"],
            col_ids_test="ids var2",
            ref_path=cfg["inputs"]["03_backup_vs_noloss"]["path"],
            values_path=values_path,
            cols_value=cols_value,
            ## Defaults:
            bu_score_gt=bu_score_gt,  # 0
            kws_scoring=kws_scoring,  #  {"value_min": 0.1}
            kws_ro=kws_ro,  # {}
            ## Optional: old for back-compatible
            mod_cfgs=mod_cfgs,  # None
            mod_names=mod_names,  # None
            ## Optional: for mainly testing
            sids=sids,  # None
            ## Optional: validation
            kws_scoring_flt=kws_scoring_flt,  # None
            table_ext=table_ext,  #'tsv'
        )
    )
    test_params(params)

    scoring_dir_path = Path(cfg["ro"]["path"]).with_suffix("")
    to_dict(
        params,
        f"{scoring_dir_path}_pms.yaml",
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
