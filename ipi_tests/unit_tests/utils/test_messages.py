"""Tests the identification of the i-PI tree written in the output headers."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import pytest

from ipi.utils.messages import get_git_info


@pytest.mark.parametrize("packed", [False, True])
def test_git_info_reads_loose_and_packed_refs(tmp_path, packed):
    """The commit of the branch is found also once Git has packed its ref."""

    commit = "a" * 40
    (tmp_path / "HEAD").write_text("ref: refs/heads/bugfix/barostats\n")
    if packed:
        (tmp_path / "packed-refs").write_text(
            "# pack-refs with: peeled fully-peeled sorted\n"
            f"{commit} refs/heads/bugfix/barostats\n"
            f"{'b' * 40} refs/remotes/origin/bugfix/barostats\n"
        )
    else:
        branch = tmp_path / "refs" / "heads" / "bugfix"
        branch.mkdir(parents=True)
        (branch / "barostats").write_text(commit + "\n")
    info = get_git_info(str(tmp_path))
    assert info["branch_name"] == "bugfix/barostats"
    assert info["last_commit"] == commit
    assert info["remote_url"] == "unknown"


def test_git_info_without_the_commit(tmp_path):
    (tmp_path / "HEAD").write_text("ref: refs/heads/main\n")
    assert get_git_info(str(tmp_path))["last_commit"] == "unknown"
