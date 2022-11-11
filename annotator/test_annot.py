import run
from datetime import datetime

arguments = {
    "infile": "../tests/test_100.xlsx",
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


def test_main():

    args = arguments.copy()
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_infile1():

    args = arguments.copy()
    args["infile"] = "../tests/testSM_100.xlsx"
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_infile2():

    args = arguments.copy()
    args["level"] = "nterm"
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_software1():

    args = arguments.copy()
    args["infile"] = "../tests/testSM_100.xlsx"
    args["software"] = "sm"
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_software2():

    args = arguments.copy()
    args["software"] = "pd"
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_sleeptime():

    args = arguments.copy()
    args["sleeptime"] = 0
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_extra():

    args = arguments.copy()
    args["noexo"] = True
    args["nomerops"] = True
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_output1():

    args = arguments.copy()
    args["outfile_type"] = "csv"
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)


def test_output2():

    args = arguments.copy()
    args["outfile_type"] = "pkl"
    args["timestamp"] = datetime.now().strftime(format="%d%m%y%H%M%S")
    run.main(args)
