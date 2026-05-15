#!/bin/bash

PACKAGES=("linalg" "optimize" "stats" "signal" "sparse" "integrate" "interpolate")

mkdir -p "CProfileResults"

for PKG in "${PACKAGES[@]}"; do
    IGNORE_ARGS=""
    if [ "$PKG" == "linalg" ]; then
        IGNORE_ARGS="--ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/tests/_cython_examples"
    elif [ "$PKG" == "integrate" ]; then
        IGNORE_ARGS="--ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/tests/test_quadrature.py"
    elif [ "$PKG" == "stats" ]; then
        IGNORE_ARGS="--ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/tests/test_continuous.py --ignore=/Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/tests/test_stats.py"
    fi

    python3 CProfileToExcel.py \
        --output "CProfileResults/scipy_${PKG}_cprofile.xlsx" \
        --package "scipy.$PKG" $IGNORE_ARGS
        
done