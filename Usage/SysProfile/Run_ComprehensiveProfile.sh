#!/bin/bash

PACKAGES=("linalg" "optimize" "stats" "signal" "sparse" "integrate" "interpolate")
OUTPUT_DIR="ComprehensiveResults"
mkdir -p $OUTPUT_DIR

for PKG in "${PACKAGES[@]}"; do
    OUTPUT_FILE="${OUTPUT_DIR}/${PKG}_Profile.xlsx"
    
    IGNORE_ARGS=""
    if [ "$PKG" == "linalg" ]; then
        IGNORE_ARGS="--ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/tests/_cython_examples"
    elif [ "$PKG" == "integrate" ]; then
        IGNORE_ARGS="--ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/tests/test_quadrature.py"
    elif [ "$PKG" == "stats" ]; then
        IGNORE_ARGS="--ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/tests/test_continuous.py --ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/tests/test_stats.py"
    fi

    python3 ComprehensiveProfile.py \
        --output "$OUTPUT_FILE" \
        --exclude-libs \
        -m pytest \
        $IGNORE_ARGS \
        --pyargs "scipy.$PKG"
done

python3 MergeComprehensiveResults.py \
    --pattern "${OUTPUT_DIR}/*_Profile.xlsx" \
    --output "${OUTPUT_DIR}/Comprehensive_Profile_Merged.xlsx"
