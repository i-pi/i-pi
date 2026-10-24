#!/bin/bash

set -e

python -c "import upet; upet.save_upet(model='pet-mad', size='xs', version='1.6.0', output='pet-mad-xs-v1.6.0.pt')"
