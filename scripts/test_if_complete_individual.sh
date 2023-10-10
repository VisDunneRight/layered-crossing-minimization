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
    count=$(tail -n 200 "$file" | cut -d ',' -f 11 | grep -oE '[0-9]+(\.[0-9]+)?' | awk '$1 > 300 {count++} END {print count}')
    if [[ "$count" -gt 50 ]]; then
        # Get the 5th comma-separated value of the 200th line from the end
        value_200th_line=$(tail -n 200 "$file" | head -n 1 | cut -d ',' -f 4 | grep -oE '[0-9]+')
        value_200th_line_truncated=$(($value_200th_line / 10 * 10))

        # Get the 5th comma-separated value of the last line
        value_last_line=$(tail -n 1 "$file" | cut -d ',' -f 4 | grep -oE '[0-9]+')
        value_last_line_truncated=$(($value_last_line / 10 * 10))
        if [ "$value_200th_line_truncated" -ne "$value_last_line_truncated" ]; then
            echo "$file : has not completed current bucket $value_last_line_truncated"
            ((incomplete_count++))
        else
            echo "$file : $count / 200 failed, bucket $value_last_line_truncated"
            ((complete_count++))
        fi
    else
        value_200th_line=$(tail -n 200 "$file" | head -n 1 | cut -d ',' -f 4 | grep -oE '[0-9]+')
        value_200th_line_truncated=$(($value_200th_line / 10 * 10 ))
        value_last_line=$(tail -n 1 "$file" | cut -d ',' -f 4 | grep -oE '[0-9]+')
        value_last_line_truncated=$(($value_last_line / 10 * 10 ))
        if [ "$value_200th_line_truncated" -ne "$value_last_line_truncated" ]; then
            echo "$file : has not completed current bucket $value_last_line_truncated"
        else
            echo "$file : $count / 200 failed, bucket $value_last_line_truncated"
        fi
        ((incomplete_count++))
    fi
done

echo "Incomplete: $incomplete_count"
echo "Complete: $complete_count"