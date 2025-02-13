# !/bin/bash

for value in {0..25}
do
	python3 benchmark_exec.py $value &
done
