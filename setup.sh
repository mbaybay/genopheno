#!/usr/bin/env bash

# Create environment:
conda create -n genopheno_py27 python=2.7 -f environment.yml

# Start environment:
source activate genopheno_py27

# Update environment
conda update