#!/usr/bin/env python
# coding: utf-8
# # Filtering
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    col_group=None,
    col_id=None,
    col_sid=None,
    cols=None,
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

    logging.basicConfig(level=logging.INFO)
    ## system functions
    from pathlib import Path

    ## data functions
    ## system functions from roux
    from roux.lib.io import read_table, to_table
    ## workflow functions from roux
    # sys.path.append('..')

    ### output set-up
    output_dir_path = Path(output_path).with_suffix("")
    Path(output_dir_path).mkdir(parents=True, exist_ok=True)
    logging.info(f"Output directory: {output_dir_path}")
    # ## to collect the stats
    # ## Inputs

    if cols is not None:
        col_var1, col_var2 = cols
    del cols
    # ### Mutation effects by exons

    df0_me = read_table(input_path)
    df0_me.head(1)
    # ### By data availlbility

    df1 = (
        df0_me.log(col_group)
        .loc[
            :,
            [
                col_group,
                col_id,
                col_sid,
                col_var1,
                col_var2,
            ],
        ]
        .rd.assert_no_dups()
        .log()
        .groupby(
            [
                col_group,
                col_sid,
            ]
        )
        .filter(
            lambda df: (df[col_id].nunique() > 1 and df[col_var2].nunique() == 2),
        )
        .log(col_group)
        .sort_values([col_group, col_sid])
    )
    df1.head(1)
    # ## Mapping backup to loss

    from roux.lib.set import list2str

    df2_ = (
        df1.query(f"`{col_var2}`==True")
        .groupby(
            [
                col_sid,
                col_group,
            ]
        )
        .agg(
            **{
                "ids var2": (
                    col_id,
                    lambda x: list2str(x, fmt="id"),
                )
            }
        )
        .reset_index()
    )
    df2_.head(1)

    df2 = df2_.log.merge(
        right=df1.query(f"`{col_var2}`==False").drop([col_var2], axis=1),
        how="inner",
        on=[col_group, col_sid],
        validate="1:m",
    )
    df2.head(1)
    # ## Output data

    to_table(
        df2,
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
