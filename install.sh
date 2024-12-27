#!/bin/bash

# Define the directory containing your script
SCRIPT_DIR="$(pwd)"/cas_finder

# Choose the appropriate file (.bash_profile or .bashrc)
TARGET_FILE="$HOME/.bashrc"
if [ ! -f "$TARGET_FILE" ]; then
    TARGET_FILE="$HOME/.bash_profile"
fi

# Add SCRIPT_DIR to PATH if not already added
if ! grep -q "export PATH=.*:$SCRIPT_DIR" "$TARGET_FILE"; then
    echo "Adding $SCRIPT_DIR to PATH in $TARGET_FILE..."
    echo "export PATH=\"\$PATH:$SCRIPT_DIR\"" >> "$TARGET_FILE"
    echo "Added successfully!"
else
    echo "$SCRIPT_DIR is already in PATH."
fi

# Apply the changes
echo "Sourcing $TARGET_FILE..."
source "$TARGET_FILE"

