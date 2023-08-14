#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

directory="$1"

if [ ! -d "$directory" ]; then
    echo "Error: $directory is not a valid directory."
    exit 1
fi

incomplete_count=0
complete_count=0

for file in "$directory"/*.csv; do
    count=$(tail -n 100 "$file" | cut -d ',' -f 11 | grep -oE '[0-9]+(\.[0-9]+)?' | awk '$1 > 60 {count++} END {print count}')
    if [ "$count" -gt 50 ]; then
        ((complete_count++))
    else
        ((incomplete_count++))
    fi
done

echo "Incomplete: $incomplete_count"
echo "Complete: $complete_count"
