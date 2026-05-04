@echo off
set "CATALYSIS_TOOLKIT_START_PATH=/xrd"
title Catalysis Data Toolkit - DEBUG MODE
echo.
echo  ============================================
echo   Catalysis Data Toolkit - DEBUG MODE
echo   (Window stays open on exit)
echo  ============================================
echo.

:: ── Try conda first (needed for GSAS-II) ────────────────────────────────
echo  [DEBUG] Step 1: Checking for conda...
where conda
if errorlevel 1 (
    echo  [DEBUG] conda NOT found - falling back to venv
    goto :USE_VENV
)
echo  [DEBUG] conda found.
echo.

:: ── CONDA PATH ──────────────────────────────────────────────────────────
echo  [Conda detected - using conda environment]
echo.

:: Get full path to this directory (handles spaces)
set "TOOLKIT_DIR=%~dp0"
set "CONDA_ENV=%TOOLKIT_DIR%.conda_env"
echo  [DEBUG] Toolkit dir: %TOOLKIT_DIR%
echo  [DEBUG] Conda env:   %CONDA_ENV%
echo.

:: Create conda environment if it doesn't exist
if not exist "%CONDA_ENV%\python.exe" (
    if not exist "%CONDA_ENV%\Scripts\python.exe" (
        echo  [DEBUG] Conda env not found - creating with Python 3.13...
        echo  Creating conda environment...
        echo  ^(First-time setup - this takes a few minutes^)
        echo.
        call conda create -p "%CONDA_ENV%" python=3.13 -y
        if errorlevel 1 (
            echo  [DEBUG] ERROR: Failed to create conda environment.
            echo  Falling back to venv...
            goto :USE_VENV
        )
        echo  [DEBUG] Conda environment created.
        echo.
    )
)

:: Activate conda environment
echo  [DEBUG] Step 2: Activating conda environment...
call conda activate "%CONDA_ENV%"
if errorlevel 1 (
    echo  [DEBUG] ERROR: Failed to activate conda environment.
    echo  Try running: conda init cmd.exe
    echo  Then close and reopen this window.
    echo  Falling back to venv...
    goto :USE_VENV
)
echo  [DEBUG] Conda env activated.
echo.

:: Show which python we're using
echo  [DEBUG] Step 3: Python location and version:
where python
python --version
echo.

:: ── Check if GSAS-II is already installed ───────────────────────────────
echo  [DEBUG] Step 4: Checking GSAS-II...
python -c "import GSASII.GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo  [DEBUG] GSAS-II found (new-style import)
    python -c "exec('try:\n import GSASII.GSASIIscriptable as G\n print(\"  GSAS-II binary dir:\", G.G2path.get_BinaryPrefix())\nexcept: print(\"  (could not determine binary dir)\")')"
    goto :GSAS_DONE
)
python -c "import GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo  [DEBUG] GSAS-II found (legacy import)
    goto :GSAS_DONE
)

:: ── GSAS-II not found — try to install ──────────────────────────────────
echo  [DEBUG] GSAS-II NOT found - attempting install...
echo  Installing GSAS-II ^(one-time, may take several minutes^)...
echo.

:: Try conda package gsas2pkg first
echo  [DEBUG] Trying: conda install gsas2pkg -c briantoby ...
call conda install gsas2pkg -c briantoby -y
if errorlevel 1 goto :GSAS_TRY_GIT

:: Install GSAS-II's optional dependencies needed for CIF import
echo  [DEBUG] Installing pycifrw and xmltodict...
call pip install pycifrw xmltodict
echo.

:: Verify conda install worked
python -c "import GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo  [DEBUG] GSAS-II installed successfully via conda.
    echo.
    goto :GSAS_DONE
)
python -c "import GSASII.GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo  [DEBUG] GSAS-II installed successfully via conda.
    echo.
    goto :GSAS_DONE
)

:GSAS_TRY_GIT
:: Conda failed — try git clone + pip install
where git >nul 2>&1
if errorlevel 1 (
    echo.
    echo  [DEBUG] WARNING: GSAS-II conda install failed and git not found.
    echo  To install GSAS-II manually, run: conda install gsas2pkg -c briantoby
    echo  Or install git from https://git-scm.com/downloads and re-run.
    echo  GSAS-II button will be hidden, but Le Bail and Rietveld still work.
    echo.
    goto :GSAS_DONE
)

echo  [DEBUG] Conda package failed, trying GitHub install...
echo.
set "GSAS_CLONE=%~dp0.gsas2_src"
if not exist "%GSAS_CLONE%" (
    git clone --depth 1 https://github.com/AdvancedPhotonSource/GSAS-II.git "%GSAS_CLONE%"
)
call pip install "%GSAS_CLONE%"

:: Install GSAS-II's optional dependencies needed for CIF import
call pip install pycifrw xmltodict
echo.

:: Verify git install worked
python -c "import GSASII.GSASIIscriptable" >nul 2>&1
if not errorlevel 1 (
    echo.
    echo  [DEBUG] GSAS-II installed successfully via pip.
    echo.
    goto :GSAS_DONE
)

echo.
echo  [DEBUG] WARNING: GSAS-II installation failed.
echo  GSAS-II button will be hidden, but Le Bail and Rietveld still work.
echo.

:GSAS_DONE
echo.

:: ── Install pip dependencies if missing ─────────────────────────────────
echo  [DEBUG] Step 5: Checking Python dependencies...
python -c "import flask, yaml, numpy, scipy, matplotlib, requests, pymatgen, pandas" >nul 2>&1
if errorlevel 1 (
    echo  [DEBUG] Some dependencies missing - installing...
    echo  Installing Python dependencies...
    echo  ^(First run: pymatgen is ~500 MB - please be patient^)
    echo.
    call pip install flask pyyaml numpy pandas scipy matplotlib requests pymatgen pycifrw xmltodict openpyxl
    if errorlevel 1 (
        echo.
        echo  [DEBUG] ERROR: Dependency installation failed.
        echo  Check your internet connection and try again.
        pause
        exit /b
    )
    echo.
    echo  [DEBUG] Dependencies installed successfully.
    echo.
) else (
    echo  [DEBUG] All dependencies present.
)
echo.

:: Verify individual deps for diagnostic output
echo  [DEBUG] Step 6: Individual dependency check:
python -c "import flask; print('  flask', flask.__version__)"
python -c "import numpy; print('  numpy', numpy.__version__)"
python -c "import scipy; print('  scipy', scipy.__version__)"
python -c "import pymatgen; print('  pymatgen', pymatgen.__version__)"
python -c "import CifFile; print('  pycifrw OK')"
python -c "import xmltodict; print('  xmltodict OK')"
python -c "import yaml; print('  pyyaml OK')"
python -c "import pandas; print('  pandas', pandas.__version__)"
echo.

goto :START_APP

:: ── VENV FALLBACK (no conda available) ──────────────────────────────────
:USE_VENV
echo  [Using Python venv - install Miniforge for GSAS-II support]
echo.

:: Check Python is available
echo  [DEBUG] Checking for Python...
python --version >nul 2>&1
if errorlevel 1 (
    echo  [DEBUG] ERROR: Python not found.
    echo  Please install Python from https://www.python.org/downloads/
    echo  Make sure to check "Add Python to PATH" during install.
    pause
    exit /b
)
where python
python --version
echo.

:: Create virtual environment if it doesn't exist
if not exist ".venv\Scripts\activate.bat" (
    echo  [DEBUG] Creating virtual environment in .venv\ ...
    echo  ^(This is a one-time setup - takes 1-2 minutes^)
    echo.
    python -m venv .venv
    if errorlevel 1 (
        echo  [DEBUG] ERROR: Failed to create virtual environment.
        pause
        exit /b
    )
    echo  [DEBUG] Virtual environment created.
    echo.
)

:: Activate the virtual environment
echo  [DEBUG] Activating venv...
call .venv\Scripts\activate.bat
echo  [DEBUG] venv activated.
echo.

:: Install / upgrade dependencies if any are missing
python -c "import flask, yaml, numpy, scipy, matplotlib, requests, pymatgen, pandas" >nul 2>&1
if errorlevel 1 (
    echo  [DEBUG] Installing dependencies into .venv\ ...
    echo  ^(First run: pymatgen is ~500 MB - please be patient^)
    echo.
    call pip install flask pyyaml numpy pandas scipy matplotlib requests pymatgen pycifrw xmltodict openpyxl
    if errorlevel 1 (
        echo.
        echo  [DEBUG] ERROR: Dependency installation failed.
        echo  Check your internet connection and try again.
        pause
        exit /b
    )
    echo.
    echo  [DEBUG] Dependencies installed successfully.
    echo.
) else (
    echo  [DEBUG] All dependencies present.
)
echo.

:: ── LAUNCH ──────────────────────────────────────────────────────────────
:START_APP
echo  [DEBUG] Step 7: Starting Catalysis Data Toolkit...
echo  Opening browser at http://localhost:5000
echo.
echo  ^(Close this window to stop the server^)
echo.

python app.py

echo.
echo  ============================================
echo  DONE - Check messages above for errors.
echo  If the server crashed, scroll up to find
echo  the Python traceback.
echo  ============================================
echo.
pause
