#!/bin/bash
module swap intel/15.0.0
icc -qopenmp -O3 -g -opt-report-phase=offload main.cpp -o memTest.x
