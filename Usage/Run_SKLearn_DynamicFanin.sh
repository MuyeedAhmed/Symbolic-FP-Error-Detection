#!/bin/bash

PACKAGES=("linear_model" "ensemble" "cluster" "decomposition" "svm" "neighbors" "tree" "metrics" "preprocessing" "model_selection")
OUTPUT_DIR="SKLearn_DynamicFanInResults"
mkdir -p $OUTPUT_DIR

for PKG in "${PACKAGES[@]}"; do
    OUTPUT_FILE="${OUTPUT_DIR}/${PKG}_FanIn.xlsx"
    echo "---------${PKG}---------"
    
    python3 DynamicFanIn.py \
        --output "$OUTPUT_FILE" \
        --exclude-libs \
        -m pytest \
        --pyargs "sklearn.$PKG"
done

python3 MergeFanInResults.py \
    --pattern "${OUTPUT_DIR}/*_FanIn.xlsx" \
    --output "${OUTPUT_DIR}/SKLearn_Dynamic_FanIn_Merged.xlsx"

