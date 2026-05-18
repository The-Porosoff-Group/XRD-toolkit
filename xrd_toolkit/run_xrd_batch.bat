@echo off
set "CATALYSIS_TOOLKIT_NO_BROWSER=1"
pushd "%~dp0.."
call run.bat --batch %*
popd
