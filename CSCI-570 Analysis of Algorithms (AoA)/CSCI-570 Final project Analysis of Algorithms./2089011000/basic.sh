#!/usr/bin/env bash
# basic.sh
g++ basic.cpp -O2 -std=c++17 -o basic
./basic "$1" "$2"
