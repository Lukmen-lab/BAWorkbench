@echo off
setlocal EnableExtensions

REM Creates or updates the local mamba/conda environment for the Workbench.
REM If mamba is not on PATH, initialize conda/mamba first or edit the call below
REM to your installation, for example:
REM call "%USERPROFILE%\miniforge3\Scripts\activate.bat"

set "ENV_NAME=ba_conjugate_workbench"
set "ENV_FILE=environment.yml"
set "TOOL="

where mamba >nul 2>nul
if not errorlevel 1 set "TOOL=mamba"

if not defined TOOL (
    where conda >nul 2>nul
    if not errorlevel 1 set "TOOL=conda"
)

if not defined TOOL (
    echo Neither mamba nor conda was found on PATH.
    echo Edit setup_env.bat to call your Miniforge/Mambaforge/Anaconda activate.bat, then rerun.
    exit /b 1
)

echo Using %TOOL% to create/update %ENV_NAME%...
call :env_exists "%TOOL%" "%ENV_NAME%"
if errorlevel 1 (
    echo Environment %ENV_NAME% does not exist. Creating it from %ENV_FILE%...
    %TOOL% env create -f "%ENV_FILE%"
) else (
    echo Environment %ENV_NAME% exists. Updating it from %ENV_FILE%...
    %TOOL% env update -n "%ENV_NAME%" -f "%ENV_FILE%" --prune
)

set "STATUS=%ERRORLEVEL%"
if not "%STATUS%"=="0" (
    echo Environment setup failed with exit code %STATUS%.
    exit /b %STATUS%
)

echo Environment %ENV_NAME% is ready.
endlocal
exit /b 0

:env_exists
%~1 env list | findstr /I /C:"%~2" >nul 2>nul
exit /b %ERRORLEVEL%
