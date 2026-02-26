#!/bin/bash
# Build MASE diagrams with Mermaid CLI
# Usage: ./build_diagrams.sh [white|transparent]

BACKGROUND="${1:-transparent}"

echo "🎨 Building MASE diagrams with background: $BACKGROUND"
echo ""

mmdc -i mase_erd.mmd -o mase_erd.png -w 800 -H 1200 -b "$BACKGROUND"
echo "  ✓ mase_erd.png"

echo ""
echo "✅ All diagrams generated!"
ls -lh mase_*.png
