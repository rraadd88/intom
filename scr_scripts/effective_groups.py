#!/usr/bin/env python
# coding: utf-8
# # Effective families
# Author: Rohan Dandage

import argh


def run(input_path=None, output_path=None, colvalue=None, idt=None, col_group=None):
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

    ## system functions
    from pathlib import Path  # noqa

    ## system functions from roux
    from roux.lib.io import read_table
    from roux.lib.io import to_table
    # sys.path.append('..')

    output_dir_path = Path(output_path).with_suffix("").as_posix()
    logging.info(f"Output directory: {output_dir_path}")
    # ## Inputs

    df3 = read_table(input_path)
    df3.head(1)

    if idt is None:
        idt = [c for c in df3 if " id " in c][0].split(" id ")[0]
    # ## Groups

    # from dups.features import merge_groups #todo
    from intom.features import merge_groups

    df4 = merge_groups(
        df3,
        colvalue=colvalue,
        idt=idt,
    )
    df4.head(1)

    if col_group != "group id":
        df4 = df4.rd.renameby_replace(
            {"group": col_group.split(" id")[0]},
        )
    df4.head(1)

    to_table(
        df4,
        f"{output_dir_path}/03_effective_groups.pqt",
    )

    df5 = (
        df4.groupby(col_group)
        .apply(
            lambda df: list(
                set(df[f"{idt} id {idt}1"].tolist() + df[f"{idt} id {idt}2"].tolist())
            )
        )
        .explode()
        .to_frame(f"{idt} id")
        .reset_index()
    )
    df5.head(1)

    to_table(
        df5,
        f"{output_dir_path}/04_effective_group_members.pqt",
    )
    # ### Group members for a pair

    def get_fam_mem(
        genes_id,
        df01,
    ):
        gene_ids = genes_id.split("--")
        df_ = df01.query(expr=f"`{idt} id` == {gene_ids}")
        assert len(df_) == 2
        assert df_[col_group].nunique() == 1
        group_id = df_[col_group].tolist()[0]
        gene_ids_fam = df01.query(
            expr=f"`{col_group}`=='{group_id}' & `{idt} id` != {gene_ids}"
        )[f"{idt} id"].tolist()
        return ";".join(gene_ids_fam)

    df6 = df4.assign(
        **{
            "group members": lambda df: df[f"{idt}s id"].apply(
                lambda x: get_fam_mem(
                    genes_id=x,
                    df01=df5,
                )
            ),
        }
    )
    df6.head(1)
    # ## Output data

    to_table(df6, output_path)
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
