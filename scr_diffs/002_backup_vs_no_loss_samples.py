#!/usr/bin/env python
# coding: utf-8
# # Backup in no-loss samples
#
#
# ```mermaid
# graph LR
# ref_ids["Per group,
# agg.ed ref by Var2"]
#     --> Samples
#
# Var1["Var1
# agg per group"]
#     --> Similar_samples
#
# Similar_samples & Samples
#     --> Mapped_samples
#
# Mapped_samples & Var2
#     --> Mapped_var2
#
# classDef dropped fill:#ffffff
# ```
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    values_path=None,
    sids=None,
    values_flt_path=None,
    col_group=None,
    col_id=None,
    col_sid=None,
    cols=None,
    col_ids_skip="ids",
    col_ids_pert="ids var2",
    ids_path=None,
    drop_ids_path=None,
    method="pearson",
    corr_min=0.9,
    samples_min=5,
    exclude_zeros=False,
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

    ## data
    import pandas as pd

    ## system functions from roux
    from roux.lib.io import read_table
    from roux.lib.io import to_table
    ## visualization functions from roux

    output_dir_path = Path(output_path).with_suffix("").as_posix()
    logging.info(f"Output directory path: {output_dir_path}")

    logging.info(cols)

    col_sid_test = f"{col_sid} {cols[1]}"

    cols_id = [
        col_sid_test,
        col_ids_skip,
    ]
    logging.info(cols_id)

    col_ssim = f"r{method[0]}"
    # ## Inputs

    df0_ids = read_table(input_path)

    if len(df0_ids) == 0:
        to_table(
            pd.DataFrame(),
            output_path,
        )
        return

    df1_ids = (
        df0_ids.loc[:, [col_sid, col_group, col_ids_skip]]
        .log.drop_duplicates()
        .log(col_ids_skip)
    )
    df1_ids.head(1)
    # ### Var1 values for the sample similarity (unfiltered)

    if values_path.endswith(".h5mu"):
        import mudata as md

        dt1 = md.read(values_path, backed=True)
        logging.info(dt1)
        mod_name = dt1.mod_names[0]
        df0_vars = dt1.mod[mod_name].to_df()
        if sids is not None:
            df0_vars = df0_vars.loc[sids, :]
        df0_vars = (
            df0_vars.melt(
                ignore_index=False,
                value_name=dt1.mod[mod_name].uns.get("col", mod_name),
            )
            .reset_index()
            .merge(
                right=dt1.mod[mod_name].var.reset_index(),  # .astype({col_group:str}),
                on=col_id,
                how="inner",
                validate="m:1",
            )
            .dropna(subset=[dt1.mod[mod_name].uns.get("col", mod_name)])
        )
        dt1.file.close()
    else:
        df0_vars = read_table(
            values_path,
            use_paths=True,
            cpus=5,
        ).astype({col_sid: "str"})
    df0_vars.head(1)

    if col_group not in df0_vars:
        assert ids_path is not None
        ## mapping missing family members
        df00 = read_table(ids_path)
        df0_vars = df00.merge(
            right=df0_vars,
            how="inner",
            on=[col_id],
            validate="1:m",
        )
    # #### Agg. by group

    df1_vars = (
        df0_vars.log(col_group)
        .log(col_sid)
        .groupby([col_group, col_sid])[cols[0]]
        .agg("sum")
        .reset_index()
        .log(col_group)
        .log(col_sid)
    )
    df1_vars.head(1)
    # ### Var1 values for mapping (filtered)

    # if not values_flt_path is None:
    df0_vars_flt = read_table(
        values_flt_path,
    )
    df0_vars_flt.head(1)
    # ## Similar samples
    # ### Matrix

    df2_vars = (
        df1_vars.log(col_group)
        .log(col_sid)
        .pivot(
            index=col_group,
            columns=col_sid,
            values=cols[0],
        )
    )
    df2_vars.head(1)

    to_table(
        df2_vars,
        f"{output_dir_path}/00_corrs_input.pqt",
    )
    del df1_vars
    # ### Correlations

    if exclude_zeros:
        import numpy as np

        df2_vars = df2_vars.replace(0, np.nan)

    # %%time
    df1_ = (
        df2_vars.corr(method=method)
        # .sum()
        ## recprocal
        .melt(
            ignore_index=False,
            value_name=col_ssim,
        )
        .rename(columns={col_sid: col_sid_test}, errors="raise")
        .reset_index()
    )
    df1_.head(1)

    _ = df1_[col_ssim].hist()

    to_table(
        df1_,
        f"{output_dir_path}/00_corrs.pqt",
    )
    # ### Filtering corr.s

    df1_flt = (
        df1_.log.query(expr=f"`{col_sid}`!=`{col_sid_test}`")
        .log.query(expr=f"`{col_ssim}`>={corr_min}")
        .sort_values([col_ssim], ascending=False)
        .log(col_sid)
        .log(col_sid_test)
    )
    df1_flt.head(1)

    ## top similar both ways
    df1 = pd.concat(
        [
            (
                df1_flt.groupby(
                    col_sid,  ## reference cell is best repreentative of these n test samples
                    sort=False,
                ).head(samples_min)
            ),
            (
                df1_flt.groupby(
                    col_sid_test,  ## test cell is best repreentative of these n reference samples
                    sort=False,
                ).head(samples_min)
            ),
        ],
        axis=0,
    ).log.drop_duplicates()
    df1.head(1)

    assert df1[col_sid].nunique() == df1[col_sid_test].nunique()
    # ## Mapping similar samples
    # ### To IDs

    df2 = (
        df1_ids.rename(
            columns={
                col_sid: col_sid_test,
            },
            errors="raise",
        )
        .log(col_sid_test)
        .log(col_ids_skip)
        .log.merge(
            right=df1,
            on=col_sid_test,
            how="inner",
        )
        .log(col_sid_test)
        ## sample id should not be one of the sample id var2
        .log(col_ids_skip)
        .groupby([col_ids_skip], as_index=False)
        .apply(lambda df: df.query(f"`sample id` != {df[col_sid_test].tolist()}"))
        .log(col_ids_skip)
        .log.head()
    )
    # ### Mapping to var2

    df3 = (
        df2.assign(
            **{
                col_id: lambda df: df[col_ids_skip].str.split(";"),
            }
        )
        .explode(col_id)
        .log(col_ids_skip)
        .log.drop_duplicates()
        .log.merge(
            right=df0_vars_flt.loc[:, [col_sid, col_id] + cols],
            how="inner",
            on=[
                col_sid,
                col_id,
            ],
            validate="m:1",
        )
        .rd.assert_dense(
            subset=cols_id
            + [
                col_sid,
                col_id,
            ],
        )
        .log(col_ids_skip)
    )
    df3.head(1)
    # ## Filtering
    # ### By var2: var2=False e.g. non-mutated

    if df3[cols[1]].any():
        df3 = (
            df3.log(col_ids_skip)
            .groupby(cols_id)
            .filter(lambda df: not df[cols[1]].any())
            .log(col_ids_skip)
        )
    df3.head(1)
    # ### Custom drop samples for specified groups

    if drop_ids_path is not None:
        df_drop = read_table(drop_ids_path)
        ## filtering
        df3 = df3.rd.filter_rows(
            df_drop,
            mode="drop",
        )
    df3.head(1)

    to_table(
        df3,
        f"{output_dir_path}/01_pre_flt.pqt",
    )
    # ### By number of similar samples

    df4 = (
        df3
        # .rename(
        #     },
        # )
        .loc[
            :,
            cols_id
            + [
                col_sid,
                col_ssim,
            ],
        ]
        .log.drop_duplicates(
            subset=cols_id
            + [
                col_sid,
            ]
        )
        .sort_values(col_ssim, ascending=False)
        .log()
        .groupby(cols_id)
        .head(samples_min)
        # .log()
        # .groupby(
        #         cols_id
        #     ).filter(
        # )
        .log()
        .rd.assert_no_dups(cols_id + [col_sid])
    )
    df4.head(1)
    # ## Mapping Values to filtered samples

    df5 = (
        df3.drop(
            [
                col_ssim,
                cols[1],
            ],
            axis=1,
        )
        .log.merge(
            right=df4,
            how="inner",
            on=cols_id
            + [
                col_sid,
            ],
            validate="m:1",
        )
        .rd.assert_dense(
            subset=cols_id
            + [
                col_sid,
                col_id,
            ],
        )
        .log(col_ids_skip)
    )
    df5.head(1)
    # ## Validation

    assert all(
        (
            df5.groupby(
                cols_id,
            )[col_sid].nunique()
        )
        <= samples_min
    )
    # ## Outputs

    to_table(df5, output_path)
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
