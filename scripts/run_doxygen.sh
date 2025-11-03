#!/usr/bin/env bash
set -euo pipefail

# Run doxygen (reads Doxyfile in repo root) then optionally open the HTML page for a function.
# Usage:
#   scripts/run_doxygen.sh [--open FUNCTION_NAME] [--output-dir DIR]
# Examples:
#   scripts/run_doxygen.sh          # builds docs/doxygen/html
#   scripts/run_doxygen.sh --open do_region_hit_boundary_interaction

OPEN_FUNCTION=""
OUTPUT_DIR="docs/doxygen"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --open)
      OPEN_FUNCTION="$2"; shift 2 ;;
    --output-dir)
      OUTPUT_DIR="$2"; shift 2 ;;
    --help|-h)
      echo "Usage: $0 [--open FUNCTION_NAME] [--output-dir DIR]"; exit 0 ;;
    *) echo "Unknown arg: $1"; exit 2 ;;
  esac
done

# Run doxygen
if ! command -v doxygen >/dev/null 2>&1; then
  echo "doxygen not found. Please install doxygen and dot (graphviz)." >&2
  exit 3
fi

echo "Running doxygen..."
doxygen Doxyfile

HTML_DIR="$OUTPUT_DIR/html"
if [[ ! -d "$HTML_DIR" ]]; then
  echo "Expected HTML output at $HTML_DIR not found." >&2
  exit 4
fi

if [[ -n "$OPEN_FUNCTION" ]]; then
  # Find the generated HTML file containing the function name
  echo "Searching for function '$OPEN_FUNCTION' in generated HTML..."
  # Search for files containing the function name in the HTML output
  FILE=$(grep -RIl "\b$OPEN_FUNCTION\b" "$HTML_DIR" | head -n 1 || true)
  if [[ -z "$FILE" ]]; then
    echo "Function '$OPEN_FUNCTION' not found in generated HTML. You can grep the HTML dir: grep -R \"$OPEN_FUNCTION\" $HTML_DIR" >&2
    exit 5
  fi
  echo "Found: $FILE"
  # If running in a desktop environment, open the file. Otherwise print the path.
  if command -v xdg-open >/dev/null 2>&1; then
    xdg-open "$FILE" || true
  else
    echo "Open the file in a browser: file://$FILE"
  fi
else
  echo "Doxygen HTML generated in: $HTML_DIR"
  echo "To open a function page: scripts/run_doxygen.sh --open <function_name>"
fi
