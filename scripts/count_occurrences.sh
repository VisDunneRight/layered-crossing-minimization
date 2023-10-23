#!/bin/bash

declare counts

for file in ./data\ storage/all\ switches/direct_transitivity/*; do
    if [ -f "$file" ]; then
        vue=$(tail -n 1 "$file" | awk -F',' '{print $4}')
        value=$(( vue / 10 * 10 ))
        if [ -n "$value" ]; then
            if [ -z "${counts[$value]}" ]; then
                counts["$value"]=1
            else
                ((counts["$value"]++))
            fi
        fi
    fi
done

for value in "${!counts[@]}"; do
    echo "$value: ${counts[$value]}"
done