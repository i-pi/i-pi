#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [[ -n "${IPI_REPO_ROOT:-}" ]]; then
    REPO_ROOT="${IPI_REPO_ROOT}"
else
    REPO_ROOT=$(cd "${SCRIPT_DIR}/../../../../.." && pwd)
fi
cd "${SCRIPT_DIR}"

if [[ -n "${PYTHON:-}" ]]; then
    PYTHON_BIN="${PYTHON}"
elif command -v python >/dev/null 2>&1; then
    PYTHON_BIN=$(command -v python)
else
    PYTHON_BIN=$(command -v python3)
fi

IPI_SCRIPT="${REPO_ROOT}/bin/i-pi"
PY_DRIVER_SCRIPT="${REPO_ROOT}/bin/i-pi-py_driver"
export PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

if [[ -d results ]]; then
    PREVIOUS_RESULTS="results.previous.$(date +%Y%m%d-%H%M%S).$$"
    mv results "${PREVIOUS_RESULTS}"
    echo "Preserved the previous results as ${PREVIOUS_RESULTS}"
fi
mkdir -p results
exec > >(tee -a results/run.log) 2>&1

SERVER_PID=""
CURRENT_STAGE="setup"
CURRENT_LOGS="results/run.log"

cleanup_server() {
    if [[ -n "${SERVER_PID}" ]] && kill -0 "${SERVER_PID}" 2>/dev/null; then
        kill "${SERVER_PID}" 2>/dev/null || true
        wait "${SERVER_PID}" 2>/dev/null || true
    fi
}

record_failure() {
    local status="$1"
    local line="$2"
    set +e
    cleanup_server
    {
        echo "FAILED"
        echo "stage: ${CURRENT_STAGE}"
        echo "exit_status: ${status}"
        echo "line: ${line}"
        echo "logs: ${CURRENT_LOGS}"
    } | tee results/FAILED >&2
    exit "${status}"
}

trap cleanup_server EXIT
trap 'record_failure "$?" "$LINENO"' ERR

require_free_socket() {
    local socket_path="$1"
    if [[ -e "${socket_path}" || -L "${socket_path}" ]]; then
        echo "Socket already exists: ${socket_path}" >&2
        echo "Stop its i-PI process, or remove it after confirming it is stale." >&2
        return 1
    fi
}

wait_for_socket() {
    local socket_path="$1"
    local attempt
    for attempt in $(seq 1 600); do
        [[ -S "${socket_path}" ]] && return 0
        if ! kill -0 "${SERVER_PID}" 2>/dev/null; then
            echo "The i-PI server exited before creating ${socket_path}." >&2
            wait "${SERVER_PID}"
            return 1
        fi
        sleep 0.1
    done
    echo "Timed out waiting for ${socket_path}." >&2
    return 1
}

run_macecalculator_socket() {
    CURRENT_STAGE="ffsocket_macecalculator"
    CURRENT_LOGS="results/ffsocket_macecalculator.server.log results/ffsocket_macecalculator.client.log"
    "${PYTHON_BIN}" "${IPI_SCRIPT}" input_ffsocket_macecalculator.xml \
        > results/ffsocket_macecalculator.server.log 2>&1 &
    SERVER_PID=$!
    wait_for_socket /tmp/ipi_mace-comparison-calculator
    "${PYTHON_BIN}" run_mace.py > results/ffsocket_macecalculator.client.log 2>&1
    wait "${SERVER_PID}"
    SERVER_PID=""
}

run_py_driver_socket() {
    CURRENT_STAGE="ffsocket_py_driver"
    CURRENT_LOGS="results/ffsocket_py_driver.server.log results/ffsocket_py_driver.client.log"
    "${PYTHON_BIN}" "${IPI_SCRIPT}" input_ffsocket_py_driver.xml \
        > results/ffsocket_py_driver.server.log 2>&1 &
    SERVER_PID=$!
    wait_for_socket /tmp/ipi_mace-comparison-py-driver
    "${PYTHON_BIN}" "${PY_DRIVER_SCRIPT}" \
        -u -a mace-comparison-py-driver -m mace \
        -o template=init.xyz,model=../mace.model,device=cpu,mace_kwargs=mace_kwargs.json \
        > results/ffsocket_py_driver.client.log 2>&1
    wait "${SERVER_PID}"
    SERVER_PID=""
}

run_ffdirect() {
    CURRENT_STAGE="ffdirect"
    CURRENT_LOGS="results/ffdirect.log"
    "${PYTHON_BIN}" "${IPI_SCRIPT}" input_ffdirect.xml > results/ffdirect.log 2>&1
}

CURRENT_STAGE="preflight"
CURRENT_LOGS="results/preflight.log"
{
    echo "Python: ${PYTHON_BIN}"
    "${PYTHON_BIN}" -c \
        'import ase, numpy, torch; from mace.calculators import MACECalculator; print("Python dependencies: OK")'
    test -f "${IPI_SCRIPT}"
    test -f "${PY_DRIVER_SCRIPT}"
} > results/preflight.log 2>&1

if [[ ! -s ../mace.model ]]; then
    CURRENT_STAGE="model_download"
    CURRENT_LOGS="results/run.log"
    echo "Downloading ../mace.model"
    (cd .. && bash getmodel.sh)
fi

require_free_socket /tmp/ipi_mace-comparison-calculator
require_free_socket /tmp/ipi_mace-comparison-py-driver

run_macecalculator_socket
run_py_driver_socket
run_ffdirect

touch results/SUCCESS
echo "All runs completed. Results are in ${SCRIPT_DIR}/results"
