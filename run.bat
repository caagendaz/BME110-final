@echo off
REM BioQuery Local - Start Script (Windows)
REM Runs the Streamlit app with proper path setup

cd /d "%~dp0"
set PYTHONPATH=%CD%\src;%PYTHONPATH%
streamlit run src/app.py %*
pause
