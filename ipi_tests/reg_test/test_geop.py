import subprocess as sp
from pathlib import Path
import pytest
import time
import os
import shutil
import tempfile
from ipi_tests.reg_test.runstools import Runner

parent = Path(__file__).parent

cmd1_cmd2_folder_output = [
    ["i-pi input.xml ", "i-pi-driver -h localhost -p 33334 -m ch4hcbe", "geop/bfgs", "min.out", ],
    ["i-pi input.xml ", "i-pi-driver -u -m ch4hcbe", "geop/sd", "min.out", ],
]


@pytest.mark.parametrize("cmd1,cmd2,folder,file", cmd1_cmd2_folder_output)
def test_cmd_and_files(cmd1, cmd2, folder, file):
    runner = Runner(parent)
    runner._run(cmd1, cmd2, folder)


if __name__ == "__main__":

    for (cmd1, cmd2, folder, file) in cmd1_cmd2_folder_output:
        print("Running {} ".format(folder))
        test_cmd_and_files(cmd1, cmd2, folder, file)
