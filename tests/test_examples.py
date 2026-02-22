def set_wd():
    import os
    print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())
    if not os.getcwd().endswith("/examples"):
        os.chdir('./examples/')
        print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())

set_wd()

from roux.lib.sys import run_com
def test_diff():
    result = run_com(
        "intom diff --input-path inputs/diff.yaml"
    )
    print(result)

def test_pair_n_diff():
    result = run_com(
        "intom pair --input-path inputs/pair.yaml"
    )
    print(result)
    # result = run_com(
    #     "intom diff --input-path outputs/pair/03_robust_pms.yaml"
    # )
    # print(result)
