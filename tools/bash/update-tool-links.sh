#!/bin/bash

# This is to be run from the root directory and will create symlinks to
# all the tools under `tools/py` in `bin`, with the `i-pi-` prefix.

for tool in `ls tools/py`
do
    ln -s ../tools/py/$tool bin/i-pi-${tool%.*}
done
