#!/usr/bin/env python
# coding: utf-8
# # Robustness summary
#
#     robustness: compensation & upregulation
#     upregulation
#     downregulation
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    col_score="LFC",
    col_group=None,
    col_id=None,
    col_sid=None,
    cols_value=None,
    col_comp="compensation for all ids",
    col_robust="robustness",
    compen_min=0.5,
    value_min=None,
    volcano_path=None,
    gene_ids_exclude=None,
    func=None,
    pval_max=None,
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
    import logging
    import sys
    from pathlib import Path
    import numpy as np
    from roux.lib.io import read_table, to_table
    from roux.viz.io import to_plot

    sys.path.append("../")

    ### to ensure all required parameters are defined above
    output_dir_path = str(Path(output_path).with_suffix("")) + "/"
    logging.info(f"Output directory: {output_dir_path}")

    cols_id = [
        col_group,
        f"ids {cols_value[1]}",
    ]
    # ## Inputs
    # ### Compensation

    df0_compen = read_table(
        input_path,
    )
    df0_compen.iloc[0, :]

    if col_comp not in df0_compen:
        col_comp = col_comp.replace(" all ", " sim ")
        logging.warning(col_comp)
        assert col_comp in df0_compen, col_comp

    if value_min is not None:
        assert all(df0_compen["mean test"] >= value_min)
    # ## Scoring

    def get_robust(
        compe,
        ratio=None,
        lfc=None,
        lfc_range=[0, 1],
        compe_range=[0, 1],
    ):
        if lfc is None:
            ## transform
            lfc = np.log2(ratio)

        ## capping
        if lfc < lfc_range[0]:
            lfc = lfc_range[0]
        elif lfc > lfc_range[1]:
            lfc = lfc_range[1]
        if compe < compe_range[0]:
            compe = compe_range[0]
        elif compe > compe_range[1]:
            compe = compe_range[1]

        ## combined
        return lfc * compe

    # for r,c in [
    #     (0,0),
    #     (0,1),
    #     (1,0),
    #     (1,1),
    #     (0,2),
    #     (2,0),
    #     (2,2),
    #     (2,0.5),
    #     ]:
    #     logging.info(
    #     )

    df1 = df0_compen.assign(
        **{
            col_robust: lambda df: df.apply(
                lambda x: get_robust(
                    compe=x[col_comp],
                    lfc=x[col_score] if "LFC" in x else None,
                    ratio=x[col_score] if "LFC" not in x else None,
                ),
                axis=1,
            ),
        }
    )
    df1.head(1)
    # ### Validations

    corrs = (
        df1.loc[:, [col_score, col_comp, col_robust]]
        .replace([0, np.inf, -np.inf], np.nan)
        # .map(lambda x: np.nan if x>10 else x)
        # .plot.scatter(
        # )
        .corr(
            method="spearman",
            min_periods=10,
        )["robustness"]
        # .melt()
        #     ['value']
    )
    if len(corrs) > 3:
        assert all(corrs > 0) or corrs.isnull().all(), corrs
    logging.info(corrs)

    df1[col_robust].notna().value_counts()

    _ = df1[col_robust].hist()

    kws_plot = dict(
        plot=dict(
            x=col_robust,
            y="Significance\n(-Log$_\\mathrm{10}$($p$))",
            alpha=0.7,
        ),
    )

    from roux.stat.transform import log_pval

    data = (
        df1.rd.check_na(
            subset=["P"],
            out=False,
        )
        .fillna({"P": 1})
        .log.dropna(
            subset=[
                kws_plot["plot"]["x"],
            ]
        )
        .assign(
            **{
                kws_plot["plot"]["y"]: lambda df: df["P"].apply(lambda x: log_pval(x)),
            }
        )
        .assign(
            **{
                kws_plot["plot"]["y"]: lambda df: df[kws_plot["plot"]["y"]].replace(
                    np.inf, df[kws_plot["plot"]["y"]].replace([np.inf], np.nan).max()
                ),
            }
        )
    )
    data.head(1)

    ax = data.plot.scatter(
        **kws_plot["plot"],
    )
    try:
        from roux.viz.annot import annot_side_curved

        _ = annot_side_curved(
            data.sort_values(
                [
                    kws_plot["plot"]["y"],
                    kws_plot["plot"]["x"],
                ],
                ascending=False,
            ).head(3),
            colx=kws_plot["plot"]["x"],
            coly=kws_plot["plot"]["y"],
            col_label=col_group,
            off=0.5,  ## data coord.
            limf=[0, 0.95],  ## data coord. fraction
            loc="right",
            test=False,
        )
    except:
        ## maybe unexpected range ..
        pass
    to_plot(f"{output_dir_path}/scatter_log_pval_vs_robustness.png")
    # ## Outputs

    to_table(
        df1,
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
