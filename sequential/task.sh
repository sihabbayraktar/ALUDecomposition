#!/bin/bash

echo "Compiling...."
make -f Makefile
echo "Compiling is done."


echo "PROCESSING...."
./lu_seq 16
./lu_seq 32
./lu_seq 64
./lu_seq 128
./lu_seq 256
./lu_seq 512
./lu_seq 1024
./lu_seq 2048
