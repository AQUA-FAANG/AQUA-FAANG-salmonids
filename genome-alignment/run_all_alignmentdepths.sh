# This script runs make_alignmentDepth_bigwig.sh from Ssal or Omyk to each other (sub-)genome individually


# helper function to chain sbatch job submissions
function sbatch_afterok {
    # Read the output from sbatch (previous command)
    local job_output
    job_output=$(cat)

    # Extract the job ID using a regular expression
    if [[ $job_output =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
        local job_id=${BASH_REMATCH[1]}
    else
        echo "Error: Could not extract job ID from sbatch output."
        return 1
    fi

    # Submit the next job with dependency
    sbatch --dependency=afterok:${job_id} "$@"
}

TMPOUTDIR=$SCRATCH/halAlignmentDepth
FINAL_OUTDIR=data/halAlignmentDepth

mkdir -p $FINAL_OUTDIR

function OneToOne_alignmentDepth {
  # <ref genome> <chrom.sizes for ref> <target genome>
  sbatch make_alignmentDepth_bigwig.sh $TMPOUTDIR/${1}_to_${3} $2 \
    $1 --step 1 --noAncestors --targetGenomes $3 |\
    sbatch_afterok merge_bigwigs.sh $TMPOUTDIR/${1}_to_${3} $FINAL_OUTDIR/${1}_to_${3}.bw $2
}

sbatch make_alignmentDepth_bigwig.sh $TMPOUTDIR/Ssal_A_Ancestors Ssal_chrom.sizes \
  Ssal_A --step 1 --targetGenomes Anc0,Anc1,Anc2 |\
  sbatch_afterok merge_bigwigs.sh $TMPOUTDIR/Ssal_A_Ancestors $FINAL_OUTDIR/Ssal_A_Ancestors.bw Ssal_chrom.sizes
  
OneToOne_alignmentDepth Ssal_A Ssal_chrom.sizes Ssal_B
OneToOne_alignmentDepth Ssal_A Ssal_chrom.sizes Omyk_A
OneToOne_alignmentDepth Ssal_A Ssal_chrom.sizes Omyk_B
OneToOne_alignmentDepth Ssal_A Ssal_chrom.sizes Eluc

OneToOne_alignmentDepth Omyk_A Omyk_chrom.sizes Omyk_B
OneToOne_alignmentDepth Omyk_A Omyk_chrom.sizes Ssal_A
OneToOne_alignmentDepth Omyk_A Omyk_chrom.sizes Ssal_B
OneToOne_alignmentDepth Omyk_A Omyk_chrom.sizes Eluc