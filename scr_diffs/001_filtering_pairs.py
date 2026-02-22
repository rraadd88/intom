#!/usr/bin/env python
# coding: utf-8
# # Filtering id pairs by the search space of interacting ids
#
# ```mermaid
# graph LR
#
# loss["Primary transcripts combined per gene\nbased on var2 and\nthe backup transcript\nwith var2 data"]
#
# pairs_tid["Pairs of transcript ids: primary maped to potential backup"]
#     --> net["Network of the linked transcript ids"]
#
# net & loss
#     --> pairs_tsid["Combined IDs for potential primary mapped to combined IDs for the backup"]
#
# pairs_tsid & loss
#     --> pairs_tsid_flt["Filtered based on backup overlaps with the combined IDs"]
#
# classDef dropped fill:#ffffff
# ```
#
# ```mermaid
# graph LR
#
# subgraph Gene
#     subgraph group1["Comparison#1"]
#
#         %% subgraph p1["Primary"]
#             p1e1[1] --- p1e2[2]:::mu --- p1e3[3] --- p1e4[4]:::mu --- p1e5[5] --- p1e6[6] --- p1e7[7] --- p1e8[8]
#             p2e1[1] --- p2e2[2]:::mu --- p2e3[3] --- p2e4[4]:::mu --- p2e5[5] --- p2e6[6] --- p2e7[7] --- p2e8[8]
#         %% end
#
#         %% subgraph b1["Mutation-skipping"]
#             b2e1[1] --- b2e2["skipped"]:::skipped --- b2e3[3] --- b2e4["skipped"]:::skipped --- b2e5[5] --- b2e6[6] --- b2e7[7] --- b2e8[8]
#         %% end
#     end
#
#     subgraph group2["Comparison#2"]
#
#         %%subgraph p2["Primary"]
#             p3e1[1] --- p3e2[2]:::mu --- p3e5[5] --- p3e6[6]
#             p4e2[2]:::mu --- p4e5[5] --- p4e6[6]
#         %%end
#
#         %%subgraph b2["Mutation-skipping"]
#             b4e1[1] --- b4e2["skipped"]:::skipped --- b4e5[5] --- b4e6[6]
#         %%end
#     end
# end
#
# classDef skipped width:0px,height:0px;
# classDef mu fill:#ffffff
# ```
# Author: Rohan Dandage

import argh


def run(
    input_path=None,
    output_path=None,
    pairs_path=None,
    col_group=None,
    col_id=None,
    col_ids_skip="ids",
    col_ids_pert="ids var2",
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
    from roux.workflow.io import check_for_exit
    ## workflow functions from roux
    # sys.path.append('..')

    ### output set-up
    output_dir_path = Path(output_path).with_suffix("")
    Path(output_dir_path).mkdir(parents=True, exist_ok=True)
    logging.info(f"Output directory: {output_dir_path}")
    # ## to collect the stats
    # ## Inputs

    idt = col_id.split(" id")[0]
    idt

    pert_type = col_ids_pert.split()[-1]
    pert_type
    # ### Variable1

    ## for e.g. column transcript id is for the backup and its corresponding values are in the next col.
    df0_in = read_table(input_path)
    # ### Pairing similar term IDs

    ## biologically expected backups e.g. those that have sequence similarity
    df0_pairs = read_table(pairs_path)
    df0_pairs.head(1)
    # ## Network of similar IDs to fetch all the similar IDs of a given ID (with var2 = True)

    import networkx as nx

    g = nx.from_pandas_edgelist(
        df0_pairs,
        source=f"{idt} id {idt}1",
        target=f"{idt} id {idt}2",
    )
    logging.info(g)
    # ## Combined IDs for items in each group

    from roux.lib.set import flatten

    df1 = (
        df0_in.log(col_ids_pert)
        ## loss
        .loc[:, [col_ids_pert]]
        .drop_duplicates()
        ## split
        ### to list
        .assign(
            **{
                f"{col_id} {pert_type}": lambda df: df[col_ids_pert].apply(
                    lambda x: x.split(";")
                ),
            }
        )
        ### to rows
        .explode(f"{col_id} {pert_type}")
        ## corresponding potential backups (mutated not yet filtered out)
        .assign(
            **{
                col_ids_skip: lambda df: df.apply(
                    lambda x: (
                        list(g.neighbors(x[f"{col_id} {pert_type}"]))
                        if x[f"{col_id} {pert_type}"] in g.nodes
                        else []
                    ),
                    axis=1,
                )
            }
        )
        .groupby(col_ids_pert)[col_ids_skip]
        .agg(flatten)
        .reset_index()
        ## corresponding backups without the mutated
        ### e.g. if all tids contain the mutation no backup will remain
        .assign(
            **{
                col_ids_skip: lambda df: df.apply(
                    lambda x: set(x[col_ids_skip]) - set(x[col_ids_pert].split(";")),
                    axis=1,
                )
            }
        )
        .log(col_ids_pert)
        .log.head(1)
    )
    # ## Mapping combined IDs to make unit of observation: `group-items` e.g. gene id with list of transcript IDs

    df2 = df0_in.log.merge(
        right=df1,
        how="inner",
        validate="m:1",
        on=col_ids_pert,
    ).log.head()

    from roux.lib.set import list2str

    if col_id in df2:
        ## valid. skip ids
        df2 = (
            df2.assign(
                overlap=lambda df: df.apply(
                    lambda x: x[col_id] in x[col_ids_skip], axis=1
                )
            )
            .log.query(expr="`overlap`==True")
            .drop(["overlap"], axis=1)
        )
        df2 = df2.assign(
            **{
                # fcol_ids_skip: lambda df: df[fcol_ids_skip].apply(
                # )
                ## for which transcripts the data is available
                col_ids_skip: lambda df: df.groupby([col_group, col_ids_pert])[
                    col_id
                ].transform(lambda x: list2str(x, fmt="id"))
            }
        )
    else:
        df2 = df2.assign(
            ids=lambda df: df[col_ids_skip].apply(lambda x: list2str(x, fmt="id")),
        )
    df2.log(col_ids_pert).head(1)

    check_for_exit(
        df2,
        output_path,
    )
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
