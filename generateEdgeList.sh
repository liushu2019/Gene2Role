#!/bin/bash

# Function to display help information
show_help() {
    echo "Usage: $0 method csv_file [threshold]"
    echo
    echo "method      : The correlation calculation method. Valid options are 'eepis' or 'spearman'."
    echo "csv_file    : Path to the CSV file containing data for correlation calculation."
    echo "threshold   : Optional. The threshold value for the correlation calculation."
    echo "              Defaults are 0.4 for 'eepis' and 2 for 'spearman'."
    echo
    echo "Example: $0 eepis data.csv 0.5"
    echo "         $0 spearman data.csv"
}



# Check if the first argument is either 'eeisp' or 'spearman'
if [ "$1" != "eeisp" ] && [ "$1" != "spearman" ]; then
    echo "Invalid method: $1"
    echo "Valid methods are 'eeisp' or 'spearman'."
    exit 1
fi

# Check if at least two arguments (method and CSV file) are passed
if [ $# -lt 2 ]; then
    echo "Usage: $0 method csv_file [threshold]"
    exit 1
fi

# Get the calculation method, CSV file path, and potential threshold
method=$1
csv_file=$2
shift 2  # Remove the first two arguments (method and CSV file)

# Set the default threshold
if [ "$method" = "eeisp" ]; then
    threshold=${1:-0.4}
    python python/eeisp.py "$csv_file" ----threCDI "$threshold" --threEEI ""
else
    threshold=${1:-2}
    python python/spearman.py "$csv_file" "$threshold"

fi

