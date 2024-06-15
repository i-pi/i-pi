ENV_BASE_DIR=$(cd "$(dirname "$(readlink -f "$script_path")")" && pwd)
export PATH=$ENV_BASE_DIR/bin:$PATH
export PYTHONPATH=$ENV_BASE_DIR:$PYTHONPATH
export IPI_ROOT=$ENV_BASE_DIR

unset ENV_BASE_DIR
