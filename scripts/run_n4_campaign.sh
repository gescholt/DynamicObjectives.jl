#!/bin/bash
# ==============================================================================
# N=4 Testing Campaign - Local Runner
# ==============================================================================
#
# Run the n=4 testing campaign locally (not on HPC).
#
# Usage:
#     ./scripts/run_n4_campaign.sh
#     ./scripts/run_n4_campaign.sh --gn 10 --max-time 1800
#     ./scripts/run_n4_campaign.sh --models "DAISY_Ex3_with_input,LV_4D_Constrained"
#
# Options:
#     --output DIR       Output directory (default: test_results/n4_campaign_local)
#     --gn N             Grid parameter (default: 8)
#     --max-time SEC     Max seconds per model (default: 3600)
#     --models LIST      Comma-separated model names
#     --threads N        Julia threads (default: auto)
#     --help             Show this help message
#
# ==============================================================================

set -e

# Default values
OUTPUT_DIR="test_results/n4_campaign_local_$(date +%Y%m%d_%H%M%S)"
GN=8
MAX_TIME=3600
MODELS=""
THREADS=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --gn)
            GN="$2"
            shift 2
            ;;
        --max-time)
            MAX_TIME="$2"
            shift 2
            ;;
        --models)
            MODELS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --help)
            head -30 "$0" | tail -25
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Find project directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

echo "=============================================="
echo "N=4 Testing Campaign - Local Runner"
echo "=============================================="
echo "Project: $PROJECT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Configuration: GN=$GN, max_time=${MAX_TIME}s"
if [ -n "$MODELS" ]; then
    echo "Models: $MODELS"
fi
echo "=============================================="
echo ""

# Set environment variables
export CAMPAIGN_OUTPUT_DIR="$OUTPUT_DIR"
export CAMPAIGN_GN="$GN"
export CAMPAIGN_MAX_TIME="$MAX_TIME"
if [ -n "$MODELS" ]; then
    export CAMPAIGN_MODELS="$MODELS"
fi

# Set Julia threads
if [ -n "$THREADS" ]; then
    export JULIA_NUM_THREADS="$THREADS"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run campaign
echo "Starting campaign..."
echo ""

if [ -n "$THREADS" ]; then
    julia --project=. --threads="$THREADS" examples/n4_hpc_campaign.jl
else
    julia --project=. examples/n4_hpc_campaign.jl
fi

echo ""
echo "=============================================="
echo "Campaign complete!"
echo "Results saved to: $OUTPUT_DIR"
echo "=============================================="
