#!/usr/bin/env python
# coding: utf-8
# # "Backup vs backup" trace back
#
# 1. To ensure that the scoring is correct and
# 2. to plot the transcript wise expression
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    values_path=None,
    col_group=None,
    col_score="LFC",
    col_sid=None,
    cols=None,
    col_id=None,
    col_ids_test="ids var2",
    value_min=None,
    samples_min=1,
    p_lt=0.1,
    hits_n=3,
    func="ttest",
    src_dir_path="./",
    kernel=None,
    force=False,
    run_tasks_post=False,
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
    import logging  # noqa

    logging.basicConfig(level=logging.INFO)
    from pathlib import Path

    ## data functions
    import numpy as np  # noqa
    import pandas as pd  # noqa

    ## system functions from roux
    from roux.lib.io import read_table, to_table

    ## visualization functions
    import matplotlib.pyplot as plt  # noqa
    import seaborn as sns

    ## visualization functions from roux
    from roux.viz.ax_ import format_ax  # noqa
    from roux.viz.io import to_plot

    ## workflow functions from roux
    from roux.workflow.task import run_tasks  # noqa

    output_dir_path = Path(output_path).with_suffix("").as_posix()
    logging.info(f"Output directory path: {output_dir_path}")

    stat_suffix = "" if func == "ttest" else " (MWU test)"
    # ## Inputs

    cols_id = [col_group, col_ids_test]
    cols_groupby = cols_id + ["ids"]
    cols_groupby

    col_sids_count_test = "samples test count"
    col_sids_count_ref = "samples ref count"

    if isinstance(value_min, dict):
        if cols[0] in value_min:
            value_min = value_min[cols[0]]
        else:
            value_min = None
    value_min
    # ### Data

    df0_bvb = read_table(input_path)
    df0_bvb.head(1)

    if len(df0_bvb) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        return

    df0_bvb.iloc[0, :]
    # ### Validations
    # #### Comparison between mean and median based scores

    if col_score is None:
        col_score = df0_bvb.filter(regex="^ratio between.*mean$").columns.tolist()[0]
    col_score

    def rename_score(c):
        return c.replace("mean (backup vs loss-backup vs no-loss) ", "\n")

    rename = {
        c: rename_score(c) + " (log2 scale)"
        for c in df0_bvb.filter(regex="^ratio between.*$").columns.tolist()
    }
    if len(rename) == 2:
        print(rename)
        kws_plot = dict(plot=dict())
        kws_plot["plot"]["x"], kws_plot["plot"]["y"] = list(rename.values())
        kws_plot

        data = df0_bvb.rename(
            columns=rename,
            errors="raise",
        ).assign(
            **{
                list(rename.values())[0]: lambda df: df[list(rename.values())[0]].apply(
                    np.log2
                ),
                list(rename.values())[1]: lambda df: df[list(rename.values())[1]].apply(
                    np.log2
                ),
            }
        )
        print(data.head(1))
        from roux.viz.scatter import plot_scatter

        ax = plot_scatter(
            data=data,
            **kws_plot["plot"],
        )
        _ = format_ax(ax)
    # #### Dependence on sample size

    ax = df0_bvb.plot.scatter(
        x=col_sids_count_test,
        y=col_sids_count_ref,
        alpha=0.5,
    )
    from roux.viz.ax_ import set_equallim

    _ = set_equallim(ax, diagonal=True)

    if len(rename) == 2:
        if len(data) > 50:
            data = (
                data
                # .rd.get_qbins(
                .rd.get_bins(
                    col_sids_count_test,
                    bins=[0, 3, 10, 50, 1000],
                ).assign(
                    **{
                        f"{col_sids_count_test} bin": lambda df: df[
                            f"{col_sids_count_test} bin"
                        ].replace({"(0.999, 3.0]": "[1, 3]"}),
                    },
                )
            )
            # data.head(1)
            # data[f'{col_sids_count_test} bin'].unique()

            if data[f"{col_sids_count_test} bin"].isnull().sum() == 0:
                # for samples_min in [1,3,5,10]:

                # if value_min is not None:

                cols_scores = data.filter(like="ratio").columns.tolist()

                for x in cols_scores:
                    fig, ax = plt.subplots(figsize=[1.5, 1.5])
                    ax = sns.boxplot(
                        # .violinplot(
                        # .pointplo(
                        data=data,
                        y=f"{col_sids_count_test} bin",
                        x=x,
                        order=data[f"{col_sids_count_test} bin"]
                        .unique()
                        .categories.tolist()[::-1],
                        hue_order=[True, False],
                        showmeans=True,
                        showfliers=False,
                        ax=ax,
                    )
                    _ = ax.axvline(0, linestyle=":")
                    ax.set(
                        ylabel="sample size",
                    )
                    _ = format_ax(
                        ax,
                        kws_legend=dict(
                            bbox_to_anchor=[1.7, 1],
                            loc=1,
                        ),
                    )
                    # to_plot(
                    #     ax,
                    # )
    # ## Filtering
    # ### By sample sizes

    col_sid

    col_sids_count_test

    df0_bvb[col_sids_count_test]

    df1_bvb = df0_bvb.rd.assert_no_dups(
        cols_groupby,
    ).log.query(expr=f"`{col_sids_count_test}`>={samples_min}")
    df1_bvb.head(1)
    # ## Statistical significance
    # ### Volcano plot

    score_type = col_score.split(") ")[-1]
    score_type

    _ = df0_bvb["P"].hist(bins=50).set(xlabel=f"P{stat_suffix}")

    data = df1_bvb.log.dropna(subset=["P"]).assign(
        **{
            "id": lambda df: range(len(df)),
        }
    )
    if "LFC" not in data:
        data = data.assign(
            **{
                "LFC": lambda df: df[col_score].apply(np.log2),
            }
        )

    data.head(1)

    if len(data) > 1:
        from roux.viz.scatter import plot_volcano

        ax = plot_volcano(
            data,
            colx="LFC",
            coly="P",
            colindex="id",
            p_lt=p_lt,
            text_increase="n",
            text_decrease="n",
        )
        _ = ax.set(ylabel="Significance\n(-Log$_\\mathrm{10}$($p$))")
        to_plot(
            ax,
            prefix=f"{output_dir_path}/volcano_",
        )
    # ### Filtering
    # #### By values: mean var1 value in test is greater than the threshold

    df2_bvb = (
        df1_bvb.query(f"`mean test` > {value_min}")
        .query(f"`P` <= {p_lt}")
        .log.query(expr=f"`{col_score}` > 1")
    )
    df2_bvb
    # ## Output data

    to_table(
        df2_bvb,
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
