@echo off
setlocal EnableExtensions

REM Launches the local Streamlit Workbench.
REM If conda/mamba activation is unavailable, edit this script to call your
REM installation's activate.bat first, for example:
REM call "%USERPROFILE%\miniforge3\Scripts\activate.bat"

set "ENV_NAME=ba_conjugate_workbench"

where mamba >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    mamba run -n "%ENV_NAME%" baworkbench
    exit /b %ERRORLEVEL%
)

where conda >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    conda run -n "%ENV_NAME%" baworkbench
    exit /b %ERRORLEVEL%
)

echo Neither mamba nor conda was found on PATH.
echo Edit run_workbench.bat to initialize your conda/mamba installation.
exit /b 1
