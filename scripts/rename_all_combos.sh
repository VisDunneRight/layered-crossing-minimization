#!/bin/bash

source_dir="/Users/connorwilson/PycharmProjects/stratisfimal-python/dbackup/all switches/direct_transitivity"
destination_dir="/Users/connorwilson/PycharmProjects/stratisfimal-python/dbackup/all switches/direct_transitivity_new"

mkdir -p "$destination_dir"

for file in "$source_dir"/*; do
    filename=$(basename "$file")
    new_filename="${filename//8/}"
    new_filename="${new_filename//6/9}"
    new_filename="${new_filename//5/6}"
    new_filename="${new_filename//4/5}"
    new_filename="${new_filename//3/4}"
    new_filename="${new_filename//2/3}"
    new_filename="${new_filename//9/2}"
    new_filename="${new_filename//7/8}"
    new_filename="${new_filename//6/7}"
    new_filename="${new_filename//5/6}"
    new_filename="${new_filename//4/5}"
    new_filename="${new_filename//3/4}"
    new_filename="${new_filename//2/3}"
    new_filename="${new_filename//1/2}"
    new_filename="${new_filename//0/1}"
    numbers=$(echo "$new_filename" | grep -oE '[0-9]+')

    echo "$numbers"

    integer_list=()

    length=${#numbers}

    for ((i = 0; i < length; i++)); do
        number="${numbers:i:1}"
        number_int=$((number))
        integer_list+=("$number_int")
    done

    echo "$integer_list"

    # Sort the integer list
    sorted_list=($(printf '%s\n' "${integer_list[@]}" | sort -n))
    joined_string=$(printf '%s' "${sorted_list[@]}")
    new_filename="exp_${joined_string}.csv"

    mv "$file" "$destination_dir/$new_filename"
done
