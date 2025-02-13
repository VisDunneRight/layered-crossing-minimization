# !/bin/bash

for value in {0..1}
do
	python3 benchmark_exec.py $value &
done
