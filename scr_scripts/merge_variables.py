#!/usr/bin/env python
# coding: utf-8
# # Scoring input data
#
# **Binarirization of sparse variables**:
#
# 1. `cols_binary={col: True categories}` (default)
#
# - the value set by`cols_binary` -> True
# - missing values -> False
# - the value *not* set by`cols_binary` -> removed
#
# Example: `{"deleteriousness": ["high"]}`
# - high deleteriousness -> True
# - missing (no) mutation -> False (if `merge_how = 'left'`)
# - moderate/low deleteriousness -> removed
#
# 2. `cols_binary={col: {True: categories, False: categories}}`
#
# - the True categories set by`cols_binary` -> True
# - the False categories set by`cols_binary` -> False
# - missing values -> False
# - the value *not* set by`cols_binary` -> removed
#
# Example: `{"deleteriousness": {"high":True,'low':False}}`
# - high deleteriousness -> True
# - low deleteriousness -> False
# - missing (no) mutation -> False (if `merge_how = 'left'`)
# - moderate deleteriousness -> removed
#
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    input_paths=None,
    sids=None,
    mod_names=None,
    paired=False,
    cols=None,
    col_group=None,
    cols_id=None,
    col_sid="sample id",
    cols_binary={},
    cols_fillna=None,
    drop_ids_path=None,
    merge_how="left",
    values_nunique_min=None,
    values_nonzero_fraction_min=None,
    values_equal_fraction_max=None,
    samples_min=None,
    kws_read_table={},
    idt=None,
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

    logging.basicConfig(level=logging.INFO)
    from pathlib import Path
    import pandas as pd
    from roux.lib.io import read_table, to_table
    import roux.lib.dfs as rd  # noqa
    # import sys
    # sys.path.append('../')

    output_dir_path = str(Path(output_path).with_suffix("").as_posix())
    logging.info(f"Output directory: {output_dir_path}")
    # ## Input data
    # ### IDs e.g. pairs

    ## group members
    df00 = read_table(input_path)

    if idt is None:
        if paired:
            if cols_id is None:
                cols_id = [c for c in df00 if " id " in c]
            idt = cols_id[0].split(" id ")[0]
            if col_group is None:
                col_group = f"{idt}s id"
            logging.info(f"cols_id={cols_id}")
            logging.info(f"col_group={col_group}")
        else:
            idt = cols_id[0].split(" id")[0]
            assert col_group is not None, col_group
        logging.info(idt)

    _ = df00.rd.assert_no_dups(subset=[col_group] + cols_id)
    # ### Modalities

    if paired:
        reciprocal = len(set(mod_names)) == 2
    else:
        reciprocal = False
    logging.info(reciprocal)

    # if reciprocal:
    if paired:
        cols_groupby = [col_group, "suffix"]
    else:
        cols_groupby = [col_group] + cols_id

    cols_context_id = [col_sid, col_group] + cols_id
    logging.info(cols_context_id)

    cols_vars = {
        "var1": cols_context_id + [cols[0]],
        "var2": cols_context_id + [cols[1]],
    }
    logging.info(cols_vars)
    # ## Variables

    if isinstance(input_paths, str) and input_paths.endswith((".h5mu")):  # ,'.h5ad')):
        import mudata as md
        # g: adopt the new default behavior to silence the warning

        dt1 = md.read(input_paths, backed=True)
        logging.info(dt1)
        dfs = {}
        for i, mod_name in enumerate(mod_names):
            dfs[f"var{i+1}"] = dt1.mod[mod_name].to_df()
            if sids is not None:
                dfs[f"var{i+1}"] = dfs[f"var{i+1}"].loc[sids, :]
            if i == 0:
                if col_group in dt1.mod[mod_name].var.reset_index():
                    df0_ids = dt1.mod[mod_name].var.reset_index()
                else:
                    df0_ids = read_table(input_path)
                df0_ids = df0_ids.loc[:, cols_id + [col_group]].astype({col_group: str})
            expr_binary = cols_binary.get(mod_name, "")
            if "`" not in expr_binary:
                dfs[f"var{i+1}"] = (
                    dfs[f"var{i+1}"]
                    .melt(
                        ignore_index=False,
                        value_name=dt1.mod[mod_name].uns.get("col", mod_name),
                    )
                    .reset_index()
                    .log.merge(
                        right=df0_ids,
                        on=cols_id,
                        how="inner",
                        validate="m:1",
                    )
                    .dropna(subset=[dt1.mod[mod_name].uns.get("col", mod_name)])
                )
            # elif "`" in expr_binary:
            else:
                ## generate from multiple columns
                assert len(dt1.mod[mod_name].layers.keys()) > 0
                from roux.lib.str import get_fills

                layer_names = get_fills(expr_binary, "`")
                dfs_layers = {}
                for layer_name in layer_names:
                    if layer_name not in dt1.mod[mod_name].layers:
                        logging.warning(f"layer_name: {layer_name} not found.")
                        continue
                    dfs_layers[layer_name] = (
                        pd.DataFrame(
                            dt1.mod[mod_name].layers[layer_name],
                            index=dt1.mod[mod_name].obs.index,
                            columns=dt1.mod[mod_name].var.index,
                        )
                        .melt(
                            ignore_index=False,
                            value_name=layer_name,
                        )
                        .reset_index()
                        .log.head()
                    )
                from roux.lib.dfs import merge_dfs

                dfs[f"var{i+1}"] = (
                    merge_dfs(list(dfs_layers.values()), how="inner")
                    .rd.query(
                        expr=expr_binary,
                        errors=None,  ## ignore missing cols
                    )
                    .assign(**{mod_name: True})
                    .loc[:, cols_id + [col_sid, mod_name]]
                    .log.head()
                    .merge(
                        right=df0_ids,
                        on=cols_id,
                        how="inner",
                        validate="m:1",
                    )
                )
                del dfs_layers
            # else:
            #     raise ValueError(mod_name)

        dt1.file.close()
    else:
        dfs = {}
        dfs["var1"] = read_table(
            input_paths["var1"],
            **kws_read_table,
        )

        if input_paths["var1"] == input_paths["var2"]:
            dfs["var2"] = dfs["var1"].copy()
            reciprocal = False
        else:
            dfs["var2"] = read_table(
                input_paths["var2"],
                **kws_read_table,
            )

    ## convert dtypes
    for k, c in zip(dfs.keys(), cols):
        if dfs[k][c].dtype in ["object"]:
            if sorted(dfs[k][c].unique()) == [
                str(s).encode("utf-8") for s in [False, None, True]
            ]:
                logging.info(f"converting {k}/{c} object data to bool ..")
                dfs[k] = (
                    dfs[k]
                    .assign(
                        **{
                            c: lambda df: df[c].map(eval),
                        }
                    )
                    .log.dropna(subset=[c])
                )
            ## strings encoded as bytes in hd5
            if type(dfs[k][c].values[0]) == bytes:  # noqa
                dfs[k][c] = dfs[k][c].astype(str).replace("None", None)
        elif dfs[k][c].dtype in ["float32"]:
            if sorted(dfs[k][c].unique()) == [0, 1]:
                logging.info(f"converting {k}/{c} float data to bool ..")
                dfs[k] = (
                    dfs[k]
                    .assign(
                        **{
                            c: lambda df: df[c].map({0: False, 1: True}),
                        }
                    )
                    .log.dropna(subset=[c])
                )

    _ = dfs["var1"].log.head()
    _ = dfs["var2"].log.head()
    # ### Small

    for k in dfs:
        if col_group not in dfs[k]:
            dfs[k] = df00.merge(
                right=dfs[k],
                how="inner",
                on=cols_id,
                validate="1:m",
            )
        dfs[k] = dfs[k].loc[:, cols_vars[k]].rd.assert_no_dups(cols_context_id)

    _ = dfs["var1"].log.head()
    _ = dfs["var2"].log.head()
    # ### Pairing

    if paired:
        dfs["var1"] = dfs["var1"].rename(
            columns={cols[0].split(f" {idt}")[0]: cols[0]}, errors="raise"
        )
    dfs["var1"].head(1)

    if paired:
        dfs["var2"] = (
            dfs["var2"]
            .replace(
                {
                    "suffix": {
                        f"{idt}1": f"{idt}2",
                        f"{idt}2": f"{idt}1",
                    }
                }
            )
            .rename(
                columns={cols[1].split(f" {idt}")[0]: cols[1]},
                errors="raise",
            )
        )
    dfs["var2"].head(1)
    # ### Filtering by samples

    sids = list(set(dfs["var1"][col_sid].values) & set(dfs["var2"][col_sid].values))
    logging.info(len(sids))

    for k in dfs:
        dfs[k] = (
            dfs[k].log(col_sid).log.query(expr=f"`{col_sid}` == {sids}").log(col_sid)
        )
    # ### Filtering by groups

    group_ids = list(
        set(df00[col_group].tolist())
        & set(dfs["var1"][col_group].values)
        & set(dfs["var2"][col_group].values)
    )
    logging.info(len(group_ids))

    for k in dfs:
        dfs[k] = (
            dfs[k]
            .log(col_group)
            .log.query(expr=f"`{col_group}` == {group_ids}")
            .log(col_group)
        )
    # ### Filtering by both samples and groups

    if drop_ids_path is not None:
        assert merge_how in ["left", "inner"], merge_how
        ## removing from var1 is sufficient
        df_drop = read_table(drop_ids_path)

        ## filtering
        dfs["var1"] = dfs["var1"].rd.filter_rows(
            df_drop,
            mode="drop",
        )
    dfs["var1"].head(1)
    # ## Merging

    # if paired:
    merge_on = cols_groupby + [col_sid]
    # else:
    logging.info(f"merging on {merge_on}")

    if cols_fillna is not None:
        if cols[1] in cols_fillna:
            ## assuming no data -> False category
            merge_how = "left"
            logging.warning(
                f"forcing merge_how={merge_how} because {cols[1]} is sparse"
            )
            logging.info(dfs["var2"][cols[1]].value_counts())
    else:
        cols_fillna = {}

    df1 = dfs["var1"].log.merge(
        right=dfs["var2"],
        on=merge_on,
        how=merge_how,
        validate="1:1",
        suffixes=[f" {idt}1", f" {idt}2"],  ## applied iif paired
    )
    df1.head(1)

    raw_merged_path = to_table(
        df1,
        f"{output_dir_path}/01_raw_merged.pqt",
    )
    # ###  Binarisation
    #
    #     var2_bin_exprs={
    #         True:  "`deleteriousness`==['high']",
    #         ## drop unassigned low, moderate
    #         False: "`deleteriousness`==None",
    #     }
    #
    #     ## assume low del as ~wt, do not drop the unassigned
    #     ### useful when high effect mu.s associate with low effect mu.s
    #     var2_bin_exprs={
    #         True:  "`deleteriousness`==['high','high_mod']",
    #         False: "`deleteriousness`==['None','low','moderate']",
    #     }
    #

    if cols_binary is not None:
        for col, expr in cols_binary.items():
            if isinstance(expr, list):
                expr = {col: expr}

            ## all the nulls -> None
            df1[col] = df1[col].where(pd.notnull(df1[col]), None)

            df1 = (
                df1.rd.assign_bool(
                    expr=expr,
                    col=col,
                    ### Imputation
                    fillna=cols_fillna.get(col),  # if None, unassigned -> na -> dropped
                    verbose=True,
                    clean=True,
                )
                ## remove the unassigned
                .log.dropna(subset=[col])
                .astype({col: bool})
            )
            assert df1[col].isin([True, False]).all(), (
                df1[col].isin([True, False]).sum()
            )
            assert (
                df1[col].nunique() == 2
            ), f"{df1[col].unique()}; perhaps cols_paired_binary is not correct."
    df1.head(1)

    df2 = df1.log.dropna(
        subset=cols,
    ).log.head()

    to_table(
        df2,
        f"{output_dir_path}/02_binarised.pqt",
    )
    # ## Filtering

    from roux.stat.io import perc_label

    for c in cols:
        print(c, perc_label(df2[c] == 0))
    # ### Reciprocal

    # %%time
    if paired and not reciprocal:
        df2 = df2.log.query(expr=f"`suffix`=='{idt}1'")
    # ### By values

    # %%time
    # %run ../modules/filtering.py
    from intom.filtering import filter_pairsby_values

    df3 = filter_pairsby_values(
        df2,
        cols=cols,
        cols_groupby=[col_group],  # cols_groupby,
        values_nunique_min=values_nunique_min,  # =2,
        values_nonzero_fraction_min=values_nonzero_fraction_min,  # =0.5,
        values_equal_fraction_max=values_equal_fraction_max,  # =0.25,
        col_count=col_group,
    )
    df3.head(1)
    # ### By sample size

    if samples_min is not None:
        df5 = (
            df3.log(col_group)
            .groupby(cols_groupby)
            .filter(lambda df: len(df) > samples_min)
            .log(col_group)
        )
    else:
        df5 = df3.copy().log(col_group)
    df5.head(1)
    # ## Outputs

    if df5[cols[1]].dtype == "bool":
        logging.info(df5[cols[1]].value_counts())

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
