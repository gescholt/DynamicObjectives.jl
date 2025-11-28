#!/bin/bash
#SBATCH --job-name=n4_campaign
#SBATCH --output=logs/n4_campaign_%j.out
#SBATCH --error=logs/n4_campaign_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=standard

# ==============================================================================
# N=4 Testing Campaign - SLURM Submission Script
# ==============================================================================
#
# This script submits the n=4 testing campaign to an HPC cluster using SLURM.
#
# Usage:
#     sbatch scripts/submit_n4_campaign.slurm
#
# Or with custom parameters:
#     sbatch --time=12:00:00 --mem=32G scripts/submit_n4_campaign.slurm
#
# Environment Variables (set before sbatch or modify below):
#     CAMPAIGN_OUTPUT_DIR - Where to save results (default: test_results/n4_campaign_$SLURM_JOB_ID)
#     CAMPAIGN_GN         - Grid parameter (default: 8)
#     CAMPAIGN_MAX_TIME   - Max seconds per model (default: 3600)
#     CAMPAIGN_MODELS     - Comma-separated model names (default: all n=4 models)
#
# ==============================================================================

# Print job information
echo "=============================================="
echo "SLURM Job Information"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start Time: $(date)"
echo "Working Directory: $(pwd)"
echo "=============================================="

# Set up environment
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
cd "$PROJECT_DIR" || exit 1

# Create logs directory if it doesn't exist
mkdir -p logs

# Set campaign configuration with defaults
export CAMPAIGN_OUTPUT_DIR="${CAMPAIGN_OUTPUT_DIR:-test_results/n4_campaign_${SLURM_JOB_ID}}"
export CAMPAIGN_GN="${CAMPAIGN_GN:-8}"
export CAMPAIGN_MAX_TIME="${CAMPAIGN_MAX_TIME:-3600}"
# CAMPAIGN_MODELS can be set externally; if not set, the Julia script uses defaults

echo ""
echo "Campaign Configuration:"
echo "  Output Directory: $CAMPAIGN_OUTPUT_DIR"
echo "  Grid Parameter (GN): $CAMPAIGN_GN"
echo "  Max Time per Model: ${CAMPAIGN_MAX_TIME}s"
echo "  Models: ${CAMPAIGN_MODELS:-all n=4 models}"
echo ""

# Load Julia module (adjust for your HPC system)
# Common module systems:
#   module load julia/1.10
#   module load julia/1.9
# Uncomment and modify as needed for your cluster:
# module purge
# module load julia/1.10

# Check Julia is available
if ! command -v julia &> /dev/null; then
    echo "ERROR: Julia not found. Please load the Julia module or ensure it's in PATH."
    echo "Try: module load julia/1.10"
    exit 1
fi

echo "Julia version:"
julia --version
echo ""

# Set Julia threading (use available CPUs)
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
echo "Julia threads: $JULIA_NUM_THREADS"
echo ""

# Create output directory
mkdir -p "$CAMPAIGN_OUTPUT_DIR"

# Copy job script to output for reproducibility
cp "$0" "$CAMPAIGN_OUTPUT_DIR/job_script.slurm" 2>/dev/null || true

# Save environment to output
env | sort > "$CAMPAIGN_OUTPUT_DIR/environment.txt"

echo "=============================================="
echo "Starting N=4 Testing Campaign"
echo "=============================================="

# Run the campaign
julia --project="$PROJECT_DIR" \
      --threads="$JULIA_NUM_THREADS" \
      "$PROJECT_DIR/examples/n4_hpc_campaign.jl"

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "Campaign Complete"
echo "=============================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "Output Directory: $CAMPAIGN_OUTPUT_DIR"
echo ""

# List output files
if [ -d "$CAMPAIGN_OUTPUT_DIR" ]; then
    echo "Generated Files:"
    ls -la "$CAMPAIGN_OUTPUT_DIR"
fi

exit $EXIT_CODE
