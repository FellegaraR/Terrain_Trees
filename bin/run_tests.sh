#!/bin/bash

cd ..
make tests
cd bin/

echo "===[TEST 1]=== INDEX GENERATION === "
./test_index_generation
echo "===[TEST 2]=== FROM SOUP TO INDEXED MESH ==="
./test_soup
echo "===[TEST 3]=== CURVATURE ==="
./test_curvature
echo "===[TEST 4]=== BATCHED TOPOLOGICAL QUERIES ==="
./test_batched_queries
echo "===[TEST 5]=== SPATIAL QUERIES ==="
./test_spatial_queries
echo "===[TEST 6]=== TERRAIN FEATURE EXTRACTION"
./test_terrain_feature_extraction

rm test_index_generation test_soup test_curvature test_batched_queries test_spatial_queries test_terrain_feature_extraction 
