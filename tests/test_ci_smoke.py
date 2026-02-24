from clipper.run import _normalize_args


def test_normalize_args_backwards_compatibility_and_defaults():
    args = {
        "infile": "dummy.xlsx",
        "infile_type": "infer",
        "software": "infer",
        "level": "all",
        "dropna": False,
        "sleeptime": 0.2,
        "noexo": False,
        "nomerops": False,
        "singlecpu": False,
        "conditionfile": None,
        "stat": False,
        "stat_pairwise": False,
        "visualize": False,
        "logo": None,
        "pseudocounts": True,
        "outfile_type": "xlsx",
        "separate": False,
    }

    normalized = _normalize_args(args)

    assert normalized["output_filetype"] == "xlsx"
    assert normalized["calcstructure"] is None
    assert normalized["threadingcores"] == "max"
    assert normalized["multipletesting"] is False
    assert normalized["volcano_foldchange"] == 1.5
