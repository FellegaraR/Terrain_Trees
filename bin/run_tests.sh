#!/bin/bash

cd ..
make test_index_generation test_soup test_curvature test_morse test_morse_simplification 
make test_batched_queries test_spatial_queries test_terrain_feature_extraction test_multivariate_morse
cd bin/

echo "===[TEST 1]=== INDEX GENERATION === "
./test_index_generation
echo "===[TEST 2]=== FROM SOUP TO INDEXED MESH ==="
./test_soup
echo "===[TEST 3]=== CURVATURE ==="
./test_curvature
echo "===[TEST 4]=== MORSE TERRAIN FEATURE EXTRACTION"
./test_morse
echo "===[TEST 5]=== MORSE TERRAIN FEATURE SIMPLIFICATION ==="
./test_morse_simplification
echo "===[TEST 6]=== BATCHED TOPOLOGICAL QUERIES ==="
./test_batched_queries
echo "===[TEST 7]=== SPATIAL QUERIES ==="
./test_spatial_queries
echo "===[TEST 8]=== TERRAIN FEATURE EXTRACTION"
./test_terrain_feature_extraction
echo "===[TEST 9]=== MULTIVARIATE MORSE TERRAIN FEATURE EXTRACTION"
./test_multivariate_morse

rm test_index_generation test_soup test_curvature test_morse test_morse_simplification 
rm test_batched_queries test_spatial_queries test_terrain_feature_extraction test_multivariate_morse
