#!/bin/bash

# Find all .cc and .h files recursively
files=$(find . -type f \( -name "*.cc" -o -name "*.h" \))

# Check if clang-format is available
if ! command -v clang-format &> /dev/null; then
  echo "clang-format not found. Please install it."
  exit 1
fi

# Format each file using clang-format with in-place modification (-i)
# and fallback style set to none (-fallback-style=none)
for file in $files; do
  echo "Formatting file $file"
  clang-format -i -style=file "$file" || echo "Failed to format: $file"
done

echo "Formatting completed."
