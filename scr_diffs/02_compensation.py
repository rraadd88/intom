#!/usr/bin/env python
# coding: utf-8
# # Compensantion
#
# The extent to which the upregulated expression of the skipping entities in the perturbed samples can compesante for the expression of the group in the reference (non/pre-peturbed) saples.
# $$
# C = \frac{\bar{X}^{post}_{skip}}{\bar{X}^{pre}}
# $$
# Where:$C$: Compensation.
# $\bar{X}^{post}_{skip}$: Upregulated expression of the skipping entities in the perturbed samples, summed across entities, then mean across the perturbed samples.
# $\bar{X}^{pre}$: Expression in the non/pre-perturbed samples, summed across (all/similar) entities, then mean across the perturbed samples.
# <!-- $m$: Total number of (similar) entities in the group.
# $n^{post}_{pert}$ and $n^{pre}_{pert}$: Total number of perturbed and reference samples.
# $n$: Total number of samples. -->
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    exp_dir_path=None,
    ids_path=None,
    sids=None,
    col_id=None,
    col_group=None,
    col_sid="sample id",
    col_ids_test="ids var2",
    cols_value=None,
    col_score="LFC",
    score_gt=None,
    cpus=1,
    force=False,
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

    # from roux.global_imports import *
    ## logging functions
    import logging

    ## data functions
    import numpy as np
    import pandas as pd

    ## system functions
    from pathlib import Path
    import sys  # noqa

    ## system functions from roux
    from roux.lib.io import read_table
    from roux.lib.io import to_table

    ## visualization functions
    import matplotlib.pyplot as plt

    ## visualization functions from roux
    from roux.viz.ax_ import format_ax
    import roux.lib.df_apply as rd  # noqa

    def get_sumby_tids_meanby_sids(df0_exp, col_sid, col_value, expr):
        return (
            df0_exp.query(
                expr=expr,
            )
            .groupby(col_sid)[col_value]
            .sum()
            .mean()
        )

    output_dir_path = Path(output_path).with_suffix("").as_posix()
    logging.info(f"Output directory path: {output_dir_path}")

    cols_id = [col_group, col_id]
    if isinstance(cols_value, list):
        col_value = cols_value[0]
        del cols_value
    # ## Inputs
    # ### Backup upregulation

    df0_bu = read_table(
        input_path,
    )

    df0_bu = df0_bu.astype(
        {
            c: str
            for c in [
                col_group,
            ]
        }
    )

    df1_bu = df0_bu.log.query(expr=f"`{col_score}` > {score_gt}").reset_index(drop=True)

    if len(df1_bu) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        exit(0)

    df1_no_bu = df0_bu.log.query(expr=f"`{col_score}` <= {score_gt}")
    df1_no_bu.head(1)

    df1_nan = df0_bu.log.query(expr=f"`{col_score}`.isnull()").log.head()
    # ### ID.s

    if ids_path is not None:
        df0_ids = read_table(ids_path).loc[:, [col_id, col_group]]
    else:
        logging.warning(
            "since all ids not provided, restricting to: 'compensation for sim ids'"
        )
    # ### Expression

    sids = list(
        set(
            df1_bu["samples test"].str.split(";").explode().tolist()
            + df1_bu["samples ref"].str.split(";").explode().tolist()
        )
    )
    logging.info(f"len(sids)={len(sids)}")

    if Path(exp_dir_path).is_file():
        from intom.omics.io_prep import read_table_prep

        df0_exp = read_table_prep(
            exp_dir_path,
            col_group=col_group,
            col_id=col_id,
            mod_name=col_value,
            sids=sids,
        )
    elif Path(exp_dir_path).is_dir():
        df0_exp = read_table(
            [f"{exp_dir_path}/{sid}.pqt" for sid in sids],
            cpus=2,
        )

    if col_group not in df0_exp and len(df1_bu) == 1:
        ## one int means one group
        df0_exp = df0_exp.assign(**{col_group: df1_bu[col_group].tolist()[0]})

    assert col_group in df0_exp, df0_exp.columns.tolist()
    if col_group not in df0_exp and ids_path is not None:
        df0_exp = df0_exp.log.merge(
            right=df0_ids,
            on=col_id,
            how="inner",
            validate="m:1",
        )
    df0_exp = (
        df0_exp.astype(
            {
                col_id: str,
                col_group: str,
            }
        )
        .log.dropna(
            subset=[col_value],
        )
        .log.head()
    )
    # ## Expression measures
    # **Tests**

    for _, x in df1_bu.sample().iterrows():
        logging.info(
            f"{x[col_group].split(';')}: {x['ids'].split(';')} - {x[col_ids_test].split(';')}"
        )
        break

    for expr in [
        f"`{col_sid}` == {x['samples ref'].split(';')} & `{col_id}` == {x['ids'].split(';')}",
        f"`{col_sid}` == {x['samples ref'].split(';')} & `{col_id}` == {x[col_ids_test].split(';')}",
        f"`{col_sid}` == {x['samples ref'].split(';')} & `{col_group}` == {x[col_group].split(';')}",
    ]:
        logging.info(
            get_sumby_tids_meanby_sids(
                df0_exp,
                col_sid=col_sid,
                col_value=col_value,
                expr=expr,
            )
        )

    agg_funcs = {
        "backup sum mean in test": lambda x: get_sumby_tids_meanby_sids(
            df0_exp,
            col_sid=col_sid,
            col_value=col_value,
            expr=f"`{col_sid}` == {x['samples test'].split(';')} & `{col_id}` == {x['ids'].split(';')}",
        ),
        "backup sum mean in ref": lambda x: get_sumby_tids_meanby_sids(
            df0_exp,
            col_sid=col_sid,
            col_value=col_value,
            expr=f"`{col_sid}` == {x['samples ref'].split(';')} & `{col_id}` == {x['ids'].split(';')}",
        ),
        "all ids sum mean in ref": lambda x: get_sumby_tids_meanby_sids(
            df0_exp,
            col_sid=col_sid,
            col_value=col_value,
            expr=(
                ## all ids
                f"`{col_sid}` == {x['samples ref'].split(';')} & `{col_group}` == {x[col_group].split(';')}"
            ),
        ),
        "sim ids sum mean in ref": lambda x: get_sumby_tids_meanby_sids(
            df0_exp,
            col_sid=col_sid,
            col_value=col_value,
            expr=(
                ## sim ids
                f"`{col_sid}` == {x['samples ref'].split(';')} & `{col_id}` == {x[col_ids_test].split(';')+x['ids'].split(';')}"
            ),
        ),
    }
    # ### Skipping in pert. ($\bar{X}^{post}_{skip}$) and ref. samples ($\bar{X}^{pre}_{skip}$)
    #
    # Note: this is a QC step because they were already calculated during differential expression analysis
    #
    # TODO:
    #
    #     QC on few, and use the mean ref and mean test directly
    #         mean test: 'backup sum mean in test',
    #         mean ref: 'backup sum mean in ref',
    #

    outp = f"{output_dir_path}/01_bu_exp_test.pqt"
    if Path(outp).exists() and not force:
        df1_bu = read_table(
            outp,
        )

    for col_agg in [
        "backup sum mean in test",
        "backup sum mean in ref",
    ]:
        if col_agg not in df1_bu:
            logging.info(col_agg)
            df1_bu = df1_bu.assign(
                **{
                    col_agg: lambda df: (
                        df
                        # .apply(
                        .rd.apply_async(
                            cpus=cpus,
                            func=agg_funcs[col_agg],
                        )
                    ),
                }
            )

    ## test
    df_ = df1_bu.rd.assert_no_na(
        [
            "backup sum mean in test",
            "backup sum mean in ref",
        ]
    ).query("`samples test count`==1 & `ids`==1")
    assert all(df_["backup sum mean in test"] == df_["mean test"])
    assert all(df_["backup sum mean in ref"] == df_["mean ref"])
    # return df1_bu

    if not Path(outp).exists():
        to_table(
            df1_bu,
            outp,
        )

    ## post
    df1_bu = df1_bu.astype(
        {
            c: str
            for c in [
                col_group,
                "ids",
                col_ids_test,
            ]
        }
    ).log.head()
    # ### All/similar entities in ref. samples ($\bar{X}^{pre}$)

    outp = f"{output_dir_path}/03_all_var1_ref.pqt"
    if Path(outp).exists() and not force:
        df1_bu = read_table(
            outp,
        )

    for col_agg in ["all ids sum mean in ref", "sim ids sum mean in ref"]:
        if col_agg == "all ids sum mean in ref" and ids_path is None:
            continue
        logging.info(col_agg)
        df1_bu = df1_bu.assign(
            **{
                col_agg: lambda df: (
                    df
                    # .apply(
                    .rd.apply_async(
                        cpus=cpus,
                        func=agg_funcs[col_agg],
                    )
                ),
            }
        ).rd.assert_no_na([col_agg])
    # #### QCs

    ## non-expressed
    if "all ids sum mean in ref" in df1_bu:
        print(
            df1_bu.query(expr="`all ids sum mean in ref`==0").shape,
            df1_bu.query(expr="`all ids sum mean in ref`<0.1").shape,
        )
    # ## Scores

    if "all ids sum mean in ref" in df1_bu:
        df1_bu = df1_bu.assign(
            **{
                "compensation for all ids": lambda df: df.apply(
                    lambda x: (
                        np.float64(x["backup sum mean in test"])
                        / x["all ids sum mean in ref"]
                    ),
                    axis=1,
                ),
                "backup fraction in ref": lambda df: df.apply(
                    lambda x: (
                        # 100*
                        # x['compensation for all tids']
                        np.float64(x["backup sum mean in ref"])
                        / x["all ids sum mean in ref"]
                    ),
                    axis=1,
                ),
            }
        )
    df1_bu = df1_bu.assign(
        **{
            "compensation for sim ids": lambda df: df.apply(
                lambda x: (
                    np.float64(x["backup sum mean in test"])
                    / x["sim ids sum mean in ref"]
                ),
                axis=1,
            ),
        },
    ).log.head()
    # ### Validations

    data = df1_bu.copy()
    kws_plot = dict(
        plot=dict(
            y="compensation for all ids"
            if "compensation for all ids" in data
            else "compensation for sim ids"
        )
    )
    data.head(1)

    # | label: fig-scatter_compen_vs_bu
    fig, ax = plt.subplots(figsize=[2, 2])
    data.plot.scatter(x=col_score, ax=ax, **kws_plot["plot"])
    ax.axhline(
        # 100,
        1,
        linestyle=":",
        color="k",
    )
    ax.set(
        xscale="log",
    )
    _ = format_ax(
        ax,
        textwrap_width=20,
    )
    # ## Outputs
    # ### Appending no-upregulation

    df5 = pd.concat(
        [
            df1_bu,
            df1_no_bu,
            df1_nan,
        ],
        axis=0,
    ).log()
    df5.tail(1)

    assert len(df5) == len(df0_bu), (len(df5), len(df0_bu))

    to_table(
        df5,
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
