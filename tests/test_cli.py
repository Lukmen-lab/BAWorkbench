import subprocess
import sys

import pytest

from ba_conjugate_workbench import launcher


def test_cli_help():
    result = subprocess.run([sys.executable, "-m", "ba_conjugate_generator.cli", "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "generate" in result.stdout


def test_workbench_console_launcher(monkeypatch):
    captured = {}

    def fake_call(command):
        captured["command"] = command
        return 0

    monkeypatch.setattr(launcher.subprocess, "call", fake_call)

    with pytest.raises(SystemExit, match="0"):
        launcher.main()

    command = captured["command"]
    assert command[:4] == [sys.executable, "-m", "streamlit", "run"]
    assert command[4].endswith("ba_conjugate_workbench\\app.py") or command[4].endswith(
        "ba_conjugate_workbench/app.py"
    )
    assert command[-2:] == ["--browser.gatherUsageStats", "false"]
