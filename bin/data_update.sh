#!/bin/bash

# Script to update asteroid data files

DATA_DIR="src/multineas/data"
MPC_URL="https://minorplanetcenter.net/Extended_Files/nea_extended.json.gz"

echo "Updating asteroid data in $DATA_DIR..."

# Ensure data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "Creating directory $DATA_DIR..."
    mkdir -p "$DATA_DIR"
fi

# Update nea_extended.json.gz
echo "Downloading nea_extended.json.gz from $MPC_URL..."
if curl -o "$DATA_DIR/nea_extended.json.gz" "$MPC_URL"; then
    echo "Successfully updated nea_extended.json.gz"
else
    echo "Failed to download nea_extended.json.gz"
    exit 1
fi

echo "Data update complete."
