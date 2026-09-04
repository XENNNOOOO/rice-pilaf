#!/bin/bash
set -e # Exit immediately if any command fails

# Check and download ClusterONE
if [ -f "cluster_one-1.0.jar" ]; then
    printf "ClusterONE already installed.\n\n"
else
    wget https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
fi

# Check and build LazyFox
if [ -x "LazyFox" ]; then
    printf "LazyFox already installed.\n\n"
else
    git clone https://github.com/TimGarrels/LazyFox lazyfoxdir
    cd lazyfoxdir
    git reset --hard d08f3c084df19bd2a1726159f181bbe3ad6f5bf4
    
    mkdir -p build && cd build
    cmake ..
    make
    
    mv LazyFox ../../LazyFox
    cd ../..
    rm -rf lazyfoxdir
    chmod +x LazyFox
    echo "LazyFox successfully built and installed."
fi