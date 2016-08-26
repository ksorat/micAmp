#!/bin/bash
module swap intel/14.0.2
icc -openmp -O3 -g -opt-report-phase=offload main.cpp -o memTest.x
