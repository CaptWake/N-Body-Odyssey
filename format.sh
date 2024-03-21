#!/bin/bash

# Find all .cc files recursively
cc_files=$(find . -name "*.[cc|h]" -not -path "./cmake-build-debug/*/*/*")

# Check if clang-format is available
if ! command -v clang-format-17 &> /dev/null; then
  echo "clang-format not found. Please install it."
  exit 1
fi

# Format each .cc file using clang-format with Google style
for file in $cc_files; do
  clang-format-17 -i "$file"
  echo "Formatted: $file"
done

echo "All .cc files formatted!"
