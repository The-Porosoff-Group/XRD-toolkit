@echo off
title XRD Fitting Toolkit
echo.
echo  ============================================
echo   XRD Fitting Toolkit
echo  ============================================
echo.
set "CATALYSIS_TOOLKIT_START_PATH=/xrd"

:: ── Try conda first (needed for GSAS-II) ────────────────────────────────
:: Check if conda is available
where conda >nul 2>&1
if errorlevel 1 goto :USE_VENV

:: ── CONDA PATH ──────────────────────────────────────────────────────────
echo  [Conda detected - using conda environment]
echo.

:: Get full path to this directory (handles spaces)
set "TOOLKIT_DIR=%~dp0"
set "CONDA_ENV=%TOOLKIT_DIR%.conda_env"

:: Create conda environment if it doesn't exist
if not exist "%CONDA_ENV%\python.exe" (
    if not exist "%CONDA_ENV%\Scripts\python.exe" (
        echo  Creating conda environment...
        echo  ^(First-time setup - this takes a few minutes^)
        echo.
        call conda create -p "%CONDA_ENV%" python=3.13 -y -q
        if errorlevel 1 (
            echo  ERROR: Failed to create conda environment.
            echo  Falling back to venv...
            goto :USE_VENV
        )
        echo  Conda environment created.
        echo.
    )
)

:: Activate conda environment
call conda activate "%CONDA_ENV%"
if errorlevel 1 (
    echo  ERROR: Failed to activate conda environment.
    echo  Try running: conda init cmd.exe
    echo  Then close and reopen this window.
    echo  Falling back to venv...
    goto :USE_VENV
)

:: ── Check if GSAS-II is already installed ───────────────────────────────
:: gsas2pkg installs into {prefix}/GSAS-II/ which isn't on sys.path by default,
:: so we add the conda prefix paths before checking imports.
python -c "import sys,os; p=sys.prefix; [sys.path.insert(0,d) for d in [os.path.join(p,'GSAS-II','GSASII'),os.path.join(p,'GSAS-II'),os.path.join(p,'GSAS-II','backcompat')] if os.path.isdir(d)]; import GSASII.GSASIIscriptable" >nul 2>&1
if not errorlevel 1 goto :GSAS_DONE
python -c "import sys,os; p=sys.prefix; [sys.path.insert(0,d) for d in [os.path.join(p,'GSAS-II','GSASII'),os.path.join(p,'GSAS-II'),os.path.join(p,'GSAS-II','backcompat')] if os.path.isdir(d)]; import GSASIIscriptable" >nul 2>&1
if not errorlevel 1 goto :GSAS_DONE

:: ── GSAS-II not found — try to install ──────────────────────────────────
echo  Installing GSAS-II ^(one-time, may take several minutes^)...
echo.

:: Try conda package gsas2pkg first (replaces old gsas2full)
echo  Trying: conda install gsas2pkg -c briantoby ...
call conda install gsas2pkg -c briantoby -y
if errorlevel 1 goto :GSAS_TRY_GIT

:: Install GSAS-II's optional dependencies needed for CIF import
call pip install pycifrw xmltodict >nul 2>&1

:: Verify conda install worked (add gsas2pkg paths for the check)
python -c "import sys,os; p=sys.prefix; [sys.path.insert(0,d) for d in [os.path.join(p,'GSAS-II','GSASII'),os.path.join(p,'GSAS-II'),os.path.join(p,'GSAS-II','backcompat')] if os.path.isdir(d)]; import GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo  GSAS-II installed successfully via conda.
    echo.
    goto :GSAS_DONE
)
python -c "import sys,os; p=sys.prefix; [sys.path.insert(0,d) for d in [os.path.join(p,'GSAS-II','GSASII'),os.path.join(p,'GSAS-II'),os.path.join(p,'GSAS-II','backcompat')] if os.path.isdir(d)]; import GSASII.GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo  GSAS-II installed successfully via conda.
    echo.
    goto :GSAS_DONE
)

:GSAS_TRY_GIT
:: Conda failed — try git clone + pip install
where git >nul 2>&1
if errorlevel 1 (
    echo.
    echo  WARNING: GSAS-II conda install failed and git not found.
    echo  To install GSAS-II manually, run: conda install gsas2pkg -c briantoby
    echo  Or install git from https://git-scm.com/downloads and re-run.
    echo  GSAS-II button will be hidden, but Le Bail and Rietveld still work.
    echo.
    goto :GSAS_DONE
)

echo  Conda package failed, trying GitHub install...
echo.
set "GSAS_CLONE=%~dp0.gsas2_src"
if not exist "%GSAS_CLONE%" (
    git clone --depth 1 https://github.com/AdvancedPhotonSource/GSAS-II.git "%GSAS_CLONE%"
)
call pip install "%GSAS_CLONE%"

:: Install GSAS-II's optional dependencies needed for CIF import
call pip install pycifrw xmltodict >nul 2>&1

:: Verify git install worked (add gsas2pkg paths for the check)
python -c "import sys,os; p=sys.prefix; [sys.path.insert(0,d) for d in [os.path.join(p,'GSAS-II','GSASII'),os.path.join(p,'GSAS-II'),os.path.join(p,'GSAS-II','backcompat')] if os.path.isdir(d)]; import GSASII.GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo.
    echo  GSAS-II installed successfully via pip.
    echo.
    goto :GSAS_DONE
)

echo.
echo  WARNING: GSAS-II installation failed.
echo  GSAS-II button will be hidden, but Le Bail and Rietveld still work.
echo.

:GSAS_DONE

:: ── Install pip dependencies if missing ─────────────────────────────────
python -c "import flask, yaml, numpy, scipy, matplotlib, requests, pymatgen, pandas" >nul 2>&1
if errorlevel 1 (
    echo  Installing Python dependencies...
    echo  ^(First run: pymatgen is ~500 MB - please be patient^)
    echo.
    call pip install flask pyyaml numpy pandas scipy matplotlib requests pymatgen pycifrw xmltodict openpyxl
    if errorlevel 1 (
        echo.
        echo  ERROR: Dependency installation failed.
        echo  Check your internet connection and try again.
        pause
        exit /b
    )
    echo.
    echo  Dependencies installed successfully.
    echo.
)

goto :START_APP

:: ── VENV FALLBACK (no conda available) ──────────────────────────────────
:USE_VENV
echo  [Using Python venv - install Miniforge for GSAS-II support]
echo.

:: Check Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo  ERROR: Python not found.
    echo  Please install Python from https://www.python.org/downloads/
    echo  Make sure to check "Add Python to PATH" during install.
    pause
    exit /b
)

:: Create virtual environment if it doesn't exist
if not exist ".venv\Scripts\activate.bat" (
    echo  Creating virtual environment in .venv\ ...
    echo  ^(This is a one-time setup - takes 1-2 minutes^)
    echo.
    python -m venv .venv
    if errorlevel 1 (
        echo  ERROR: Failed to create virtual environment.
        pause
        exit /b
    )
    echo  Virtual environment created.
    echo.
)

:: Activate the virtual environment
call .venv\Scripts\activate.bat

:: Install / upgrade dependencies if any are missing
python -c "import flask, yaml, numpy, scipy, matplotlib, requests, pymatgen, pandas" >nul 2>&1
if errorlevel 1 (
    echo  Installing dependencies into .venv\
    echo  ^(First run: pymatgen is ~500 MB - please be patient^)
    echo.
    call pip install flask pyyaml numpy pandas scipy matplotlib requests pymatgen pycifrw xmltodict openpyxl --quiet
    if errorlevel 1 (
        echo.
        echo  ERROR: Dependency installation failed.
        echo  Check your internet connection and try again.
        pause
        exit /b
    )
    echo.
    echo  Dependencies installed successfully.
    echo.
)

:: ── LAUNCH ──────────────────────────────────────────────────────────────
:START_APP
if /I "%~1"=="--batch" goto :RUN_BATCH
if /I "%~1"=="--fetch-cifs" goto :RUN_FETCH_CIFS
echo  Starting XRD Fitting Toolkit...
echo  Opening browser at http://localhost:5000/xrd
echo.
echo  ^(Close this window to stop the server^)
echo.

python app.py

goto :APP_DONE

:RUN_BATCH
set "BATCH_ARGS="
:BUILD_BATCH_ARGS
shift /1
if "%~1"=="" goto :DO_BATCH
set BATCH_ARGS=%BATCH_ARGS% "%~1"
goto :BUILD_BATCH_ARGS
:DO_BATCH
echo  Running XRD batch workflow...
echo.
python scripts\xrd_batch.py %BATCH_ARGS%
goto :APP_DONE

:RUN_FETCH_CIFS
set "FETCH_ARGS="
:BUILD_FETCH_ARGS
shift /1
if "%~1"=="" goto :DO_FETCH
set FETCH_ARGS=%FETCH_ARGS% "%~1"
goto :BUILD_FETCH_ARGS
:DO_FETCH
echo  Fetching XRD CIF files...
echo.
python scripts\fetch_cifs.py %FETCH_ARGS%

:APP_DONE
echo.
echo  ============================================
echo  Run finished. If this was unexpected,
echo  check the error message above.
echo  ============================================
pause
