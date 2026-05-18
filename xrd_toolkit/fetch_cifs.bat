@echo off
pushd "%~dp0.."
call run.bat --fetch-cifs %*
popd
