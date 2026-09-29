#!/bin/bash

i-pi input.xml > log.i-pi &
ipi_pid=$!

wait "${ipi_pid}"
