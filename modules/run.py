import os
from pathlib import Path
import logging
# root logger level for the imported functions
# import logging
# root_logger = logging.getLogger()
# root_logger.setLevel(logging.INFO)
## for the notebook
from roux.lib.log import Logger
logging = Logger(
    # level=log_level
)


def get_src_path(
    pkg_name='intom',
    ) -> str:
    """Get the beditor source directory path.

    Returns:
        str: path
    """
    import site

    src_path = f"{site.getsitepackages()[0]}/{pkg_name}/"
    if not Path(src_path).exists():
        logging.warning(f"package is installed in development mode at {src_path}")
        # import beditor
        # src_path=str(Path(beditor.__file__).parent)
        src_path = str(Path(os.path.realpath(__file__)).parent.parent)
    else:
        logging.info(src_path)
    return src_path
    
def info():
    """
    Setup info.
    """
    get_src_path()
    

docs_common="""
Input:
    Parameters could be provided using `options` below, or a dictionary file: `--input-path params.yaml`
    Use json-string format for the non-string parameters.
    See the `{repo}/examples/inputs` for examples.
"""        
from intom.scr_diffs.run_pairing import run as pair
pair.__name__ = "pair"
pair.__doc__ = f"""(Optional) Pre-processing to pair IDs.

{docs_common}
Output: {{output_dir_path}}/03_robust.pqt

"""

from intom.scr_diffs.run_scoring import run as diff
diff.__name__ = "diff"
diff.__doc__ = f"""Obtaining robustness probabilities.

{docs_common}
Output: {{output_dir_path}}/03_robust_pms.yaml

"""

import argh

parser = argh.ArghParser()
parser.add_commands(
    [
        pair,         ## paired (ref) sids
        diff,         ## paired (ref) sids known
        # pair_n_quant, ## paired (ref) sids unknown, so inferred first
        # re_quant,     ## merged data is available
        info,
    ]
)

if __name__ == "__main__":
    
    parser.dispatch()