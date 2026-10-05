import subprocess
import sys


def test_cli_help():
    result = subprocess.run([sys.executable, "-m", "ba_conjugate_generator.cli", "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "generate" in result.stdout
