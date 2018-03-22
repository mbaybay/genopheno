#!/usr/bin/env bash

# Create environment:
conda create -n genopheno_py27 python=2.7

# Start environment:
source activate genopheno_py27

# install requirements
pip install requirements.txt

# Update environment
conda update