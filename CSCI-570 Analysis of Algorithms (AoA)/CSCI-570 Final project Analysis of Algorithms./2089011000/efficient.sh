#!/usr/bin/env bash
# efficient.sh
g++ efficient.cpp -O2 -std=c++17 -o efficient
./efficient "$1" "$2"
