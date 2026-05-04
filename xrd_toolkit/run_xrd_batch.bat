@echo off
set "CATALYSIS_TOOLKIT_NO_BROWSER=1"
call "%~dp0..\run.bat" --batch %*
