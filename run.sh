#!/bin/bash
# BioQuery Local - Start Script
# Runs the Streamlit app with proper path setup

cd "$(dirname "$0")"
export PYTHONPATH="${PWD}/src:${PYTHONPATH}"
streamlit run src/app.py "$@"
