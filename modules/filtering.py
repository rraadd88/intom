import logging
import roux.lib.dfs as rd #noqa

def filter_pairsby_values(
    df1,
    cols,
    cols_groupby,
    col_count,
    values_nunique_min=None,#=2,
    values_nonzero_fraction_min=None,#=0.5,
    values_equal_fraction_max=None,#=0.25,  
    ):
    if values_nunique_min is None:
        values_nunique_min=1
    if values_nonzero_fraction_min is None:
        values_nonzero_fraction_min=0
    if values_equal_fraction_max is None:
        values_equal_fraction_max=1
        
    for i,c in enumerate(cols):
        if df1[c].dtype=='bool':
            assert i==1, "bool values should be var2"
            logging.warning(f"{c} is boolean; values_nonzero_fraction_min disabled and values_nunique_min = 2 .")
            values_nonzero_fraction_min_=0
            values_nunique_min_=2
        else:
            values_nonzero_fraction_min_=values_nonzero_fraction_min
            values_nunique_min_=values_nunique_min
        df1=(
            df1
            .assign(
            **{
                f'var{i+1} nunique':lambda df: df.groupby(cols_groupby)[c].transform(lambda x: x.nunique()),

                f'var{i+1} nonzero':lambda df: df.groupby(cols_groupby)[c].transform(lambda x: sum(x!=0)/len(x)),

            }
            )
            # filter
            # .log
            .query(expr=f"`var{i+1} nunique`>={values_nunique_min_}")
            # clean
            .drop([f'var{i+1} nunique'],axis=1)
            # count
            .log(
                col_count,
                label=f"`var{i+1} nunique`>={values_nunique_min_}"
            )

            # filter
            # .log
            .query(expr=f"`var{i+1} nonzero`>={values_nonzero_fraction_min_}")
            # clean
            .drop([f'var{i+1} nonzero'],axis=1)
            # count
            .log(
                col_count,
                label=f"`var{i+1} nonzero`>={values_nonzero_fraction_min_}"
            )    
        )
        
    if cols!=2:
        return df1
    else:
        ## paired filtering
        return (
            df1
            .assign(
            **{
                'var1=var2':lambda df: df[cols[0]]==df[cols[1]],
                'var1=var2 fraction':lambda df: df.groupby(cols_groupby)['var1=var2'].transform(lambda x: sum(x)/len(x)),
            }
            )

        # .log
            .query(expr=f"`var1=var2 fraction`<={values_equal_fraction_max}")
            .drop(['var1=var2','var1=var2 fraction'],axis=1)
            .log(
                col_count,
                label=f"`var1=var2 fraction`<={values_equal_fraction_max}"
            )
        )