#!/usr/bin/env bash
set -euo pipefail

cxx="${CXX:-c++}"
"$cxx" -O3 -std=c++17 generate_seed_histogram.cpp -o generate_seed_histogram
./generate_seed_histogram --self-test
./generate_seed_histogram > seed_histogram_13_prime.txt
python3 verify_c45_romanoff_upper_bound.py --self-test --json-out certificate_summary.json seed_histogram_13_prime.txt
