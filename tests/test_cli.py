import pytest

from chemonym import __version__
from chemonym.__main__ import main


def test_single_smiles_prints_name(capsys):
    exit_code = main(["CC"])
    captured = capsys.readouterr()
    assert exit_code == 0
    assert captured.out == "ethane\n"


def test_multiple_smiles_prints_tab_separated(capsys):
    exit_code = main(["C", "CC"])
    captured = capsys.readouterr()
    assert exit_code == 0
    assert captured.out == "C\tmethane\nCC\tethane\n"


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])
    captured = capsys.readouterr()
    assert exc_info.value.code == 0
    assert captured.out.strip() == __version__


def test_invalid_smiles_reports_error_and_nonzero_exit(capsys):
    exit_code = main(["not_a_smiles"])
    captured = capsys.readouterr()
    assert exit_code != 0
    assert captured.out == ""
    assert "not_a_smiles" in captured.err


def test_unsupported_structure_reports_error_and_nonzero_exit(capsys):
    exit_code = main(["CCO"])
    captured = capsys.readouterr()
    assert exit_code != 0
    assert captured.out == ""
    assert "CCO" in captured.err


def test_continues_after_error(capsys):
    exit_code = main(["not_a_smiles", "CC"])
    captured = capsys.readouterr()
    assert exit_code != 0
    assert "CC\tethane" in captured.out
    assert "not_a_smiles" in captured.err
