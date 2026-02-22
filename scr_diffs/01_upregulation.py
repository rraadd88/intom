#!/usr/bin/env python
# coding: utf-8
# # Differential expression ("Backup vs backup")
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    ref_path=None,
    values_path=None,
    col_group=None,
    col_id=None,
    col_sid=None,
    cols=None,
    col_ids_test="ids var2",
    col_ids_skip_in="ids",
    col_score="LFC",
    pval_max=0.1,
    value_min=None,
    func="limma",
    col_diffid=None,
    kws_diff={"robust": True, "trend": False, "env_name": "intom"},
    exclude_zeros=True,
    include_diff_by_agg=False,
    kws_pre_pval={"plog_for_pval": True},
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
    import numpy as np
    import pandas as pd

    ## system functions from roux
    from roux.lib.io import read_table, to_table

    ### output set-up
    output_dir_path = Path(output_path).with_suffix("")
    Path(output_dir_path).mkdir(parents=True, exist_ok=True)
    logging.info(f"Output directory: {output_dir_path}")

    stat_suffix = ""  # if func in ['ttest','limma'] else " (MWU test)"

    if func == "limma":
        include_diff_by_agg = True
        assert values_path is not None
    # ## Inputs

    cols_id = [col_group, col_ids_test]
    cols_groupby = cols_id + [col_ids_skip_in]
    cols_groupby
    # ### Test samples, with pert.

    df0_bl = read_table(input_path).sort_index(axis=1)

    if "comparison" in df0_bl:
        ## direct combined input with known similar samples
        ## e.g.:  comparison	gene id	ids	ids var2	mRNA expression (TPM)	sample id	transcript id
        assert ref_path is None, ref_path
        df1 = df0_bl.rd.assert_dense(subset=df0_bl.columns.drop([cols[0]]).tolist())
    # ### Ref. samples, without pert.

    if "df1" not in locals():
        df0_bn = read_table(ref_path).sort_index(axis=1)

    if "df1" not in locals():
        ## filtering by the loss cases
        df1_bn = (
            df0_bn.log.merge(
                right=df0_bl.loc[:, cols_groupby].drop_duplicates(),
                how="inner",
                on=[col_group, col_ids_skip_in],
            )
            .loc[:, cols_groupby + [col_sid, col_id, cols[0]]]
            .log.drop_duplicates()
            .log.head()
        )
    # ## Merging

    if "df1" not in locals():
        df1 = (
            pd.concat(
                {
                    "test": df0_bl,
                    "ref": df1_bn,
                },
                axis=0,
                names=["comparison"],
            )
            .reset_index(0)
            .log()
        )
        del df0_bl, df0_bn

    ## qc
    df1 = (
        df1.groupby(cols_id)
        .filter(lambda df: df["comparison"].nunique() == 2)
        .log()
        ## validations
        .rd.assert_dense(["comparison"] + cols_id + [col_sid, col_id])
        .sort_values(cols_id)
        .log.head()
    )

    if len(df1) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        return

    to_table(df1, f"{output_dir_path}/00_combined_diff_input.pqt")
    # ### Filtering (for faster processing)

    df1_ = df1.copy()
    # #### Dropping if zeros for test and ref

    df1_ = (
        df1_.log()
        .groupby(cols_groupby)
        .filter(lambda df: df[cols[0]].sum() > 0)
        .log(col_ids_skip_in)
    )

    if len(df1_) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        return
    # #### By IDs: same number of ids for test and ref

    df1_ = (
        df1_.log()
        .groupby(cols_id)
        .filter(
            lambda df: df.groupby("comparison")[col_id].nunique().nunique() == 1,
        )
        .log(col_ids_skip_in)
    )
    df1_.head(1)

    if len(df1_) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        return
    # #### By IDs: sample ids for test and ref

    df1_ = (
        df1_.log()
        .groupby(cols_id)
        .filter(
            lambda df: df[col_ids_skip_in].nunique() == 1,
        )
        .log(col_ids_skip_in)
    )
    df1_.head(1)

    if len(df1_) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        return

    to_table(df1_, f"{output_dir_path}/01_filtered_diff_input.pqt")
    # ## Summed values by samples

    if include_diff_by_agg:
        df1_summed = (
            df1_.groupby(
                cols_groupby
                + [
                    "comparison",
                    col_sid,
                ]
            )[cols[0]]
            .sum()
            .reset_index()
            .log(col_group)
            .groupby(cols_groupby)
            .filter(
                ## not 1 sample each for test and ref
                lambda df: len(df) > 2
            )
            .log(col_group)
        )
        if (
            len(
                df1_summed.rd.check_dups(
                    subset=cols_groupby + [col_sid, cols[0]],
                )
            )
            > 0
        ):
            logging.warning(
                "dups found where same id is ref and test with the same value!"
            )
            df1_summed = (
                df1_summed.log.drop_duplicates(
                    subset=cols_groupby + [col_sid, cols[0]], keep=False
                )
                .log(col_group)
                .groupby(cols_groupby)
                .filter(
                    ## not 1 sample each for test and ref
                    lambda df: df["comparison"].nunique() == 2
                )
                .log(col_group)
            )
        print(df1_summed.head(1))
        if len(df1_summed) == 0:
            to_table(
                pd.DataFrame(),
                output_path,
            )
            return
        to_table(
            df1_summed, f"{output_dir_path}/01_filtered_diff_input_summedby_sample.pqt"
        )
    # ## Diff.

    if func == "limma":
        #     values_path
        # )
        from intom.omics.io_prep import read_table_prep

        df0_vars = read_table_prep(
            values_path,
            col_id=col_id,
            col_group=col_group,
            mod_name=cols[0],
        )

        from intom.stats.diffs import get_diff_ins, get_diff
        import roux.lib.df_apply as rd  # noqa

        # from tqdm import tqdm
        if col_diffid is None:
            logging.warning(f"col_diffid=col_group ({col_group}) ..")
        dfs_full = {}
        for k, df in df1_summed.groupby(
            by=cols_groupby,
        ):
            df0_mat, df0_sids = get_diff_ins(
                df,
                df0_vars,
                col_value=cols[0],
                col_sid=col_sid,
                col_group=col_group,
                col_diffid=col_diffid,  # if col_diffid is not: sum by col_group,
                exclude_zeros=exclude_zeros,
            )
            dfs_full[k] = get_diff(
                df0_mat,
                df0_sids,
                groups=["test", "ref"],
                col_sample=col_sid,
                col_condition="comparison",
                col_gene="gene id" if "gene id" in df0_mat else col_diffid,
                **kws_diff,
            )

        df2_full = (
            pd.concat(
                dfs_full,
                axis=0,
                names=cols_groupby,
            )
            .rename(
                columns={
                    (
                        col_group if col_diffid is None else col_diffid
                    ): f"{(col_group if col_diffid is None else col_ids_skip_in)} compared",
                    "logFC": col_score,
                    "P.Value": "P",
                },
                errors="raise",
            )
            .reset_index(list(range(len(cols_groupby))))
        )
        to_table(df2_full, f"{output_dir_path}/02_diff_output_full.pqt")
        ## filter to remove other groups compared
        df2 = df2_full.log.query(
            expr=f"`{col_group if col_diffid is None else col_ids_skip_in}` == `{col_group if col_diffid is None else col_ids_skip_in} compared`"
        )
        if len(df2) == 0:
            to_table(
                pd.DataFrame(),
                output_path,
            )
            return

        df2_means = (
            df1_summed.pivot_table(
                index=cols_groupby,
                columns="comparison",
                values=cols[0],
                aggfunc="mean",
            )
            .add_prefix("mean ")
            .reset_index()
        )

        from roux.stat.transform import get_q

        df2 = (
            df2.drop(
                [f"{col_group if col_diffid is None else col_ids_skip_in} compared"],
                axis=1,
            )
            .assign(
                change=lambda df: df[col_score].apply(
                    lambda x: "increase" if x > 0 else "decrease" if x < 0 else None
                ),
            )
            .pipe(
                get_q,
                col="P",
            )
            .log(col_group)
            .log.merge(
                right=df2_means,
                on=cols_groupby,
                how="inner",
                validate="1:1",
            )
        )
    else:
        from intom.stats.diffs import ttest

        dfs_in = {
            "not_agg": df1_,
        }
        if include_diff_by_agg:
            dfs_in["agg"] = df1_summed

        from roux.stat.diff import get_stats_groupby

        # %run ../../code/roux/roux/stat/diff.py
        dfs_out = {}
        for k, df_ in dfs_in.items():
            print(k)
            dfs_out[k] = get_stats_groupby(
                df_,
                colsubset="comparison",
                subsets=[
                    "test",
                    "ref",
                ],
                cols_value=[
                    cols[0],
                ],
                colindex=[
                    col_sid,
                    # col_ids_test
                ]
                + ([col_id] if col_id in df_ else []),
                cols_group=cols_groupby,
                func=ttest if func == "ttest" else None,  ## function
                change_type=["diff"],
                coff_p=pval_max,
                coff_q=pval_max,
            )

        df2 = pd.concat(dfs_out, axis=0, names=["agg"]).reset_index(0)

        print(
            df2["change"].value_counts(),
            df2.query("`P` < 0.1").groupby(["agg", "change"])[col_group].nunique(),
        )

        #### Ratio

        df2 = df2.assign(
            **{
                "ratio between mean (subset1-subset2)": lambda df: df["mean subset1"]
                / df["mean subset2"],
                col_score: lambda df: df["ratio between mean (subset1-subset2)"].apply(
                    lambda x: np.log2(x + 1e-3)
                ),
            }
        ).rename(
            columns={
                "mean subset1": "mean test",
                "mean subset2": "mean ref",
            },
            errors="raise",
        )

    if len(df2) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        return

    to_table(df2, f"{output_dir_path}/02_diff_output.pqt")
    # ### Stats for validation
    # #### Overall increase vs decrease

    (
        df2
        # .log.query(
        # )
        .query(expr=f"`Q{stat_suffix}`<0.1")
        .groupby("change")[col_group]
        .nunique()
    )
    # #### Common sample sizes
    # ## Post-processing
    # #### Counts of the transcripts

    df_ids = (
        df1_.loc[:, cols_id + [col_ids_skip_in]]
        .log.drop_duplicates()
        .rd.assert_no_dups(
            subset=[
                col_ids_test,
            ]
        )
        .astype(
            {
                col_ids_test: str,
            }
        )
        .assign(
            **{
                f"{col_ids_test} count": lambda df: df[col_ids_test]
                .str.split(";")
                .apply(len),
                "ids count": lambda df: df[col_ids_skip_in].str.split(";").apply(len),
            }
        )
    )
    df_ids.head(1)
    # ### Sample lists

    from roux.lib.set import list2str

    df_sids = (
        df1_.astype({col_sid: str})
        .groupby(["comparison"] + cols_id)[col_sid]
        .agg(lambda x: list2str(x, fmt="id"))
        .unstack(0)
        .add_prefix("samples ")
        .reset_index()
        .assign(
            **{
                "samples ref count": lambda df: df["samples ref"]
                .str.split(";")
                .apply(len),
                "samples test count": lambda df: df["samples test"]
                .str.split(";")
                .apply(len),
            },
        )
    )
    df_sids.head(1)

    df3 = (
        df2.astype({k: str for k in cols_groupby})
        .log.merge(
            right=df_ids.astype({k: str for k in cols_groupby}),
            on=cols_groupby,
            how="inner",
            validate="m:1",
            validate_equal_length=True,
        )
        .astype({k: str for k in cols_id})
        .log.merge(
            right=df_sids.astype({k: str for k in cols_id}),
            on=cols_id,
            how="inner",
            validate="m:1",
            validate_equal_length=True,
        )
    )
    df3.head(1)
    # ### Setting very low values to 0

    if value_min is not None:
        if "AveExpr" not in df3:
            # non-'limma'
            df3["AveExpr"] = df3.loc[:, ["mean test", "mean ref"]].mean(axis=1)
        df3 = (
            df3.rd.check_na(
                col_score,
                out=False,
            )
            .assign(
                **{
                    col_score: lambda df: df.apply(
                        lambda x: x[col_score] if x["AveExpr"] >= value_min else 0,
                        axis=1,
                    )
                }
            )
            .rd.check_na(
                col_score,
                out=False,
            )
        )
    # ## Output data

    to_table(
        df3,
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
