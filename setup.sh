#!/bin/bash

git clone https://github.com/vik16nathan/SIR_mouse_regional_tuned
cd algorithm
git clone https://github.com/srahayel/SIR_mouse/
cd ..
git clone https://github.com/yohanyee/ggslicer/ 
uv pip install --no-build-isolation -r pyproject.toml