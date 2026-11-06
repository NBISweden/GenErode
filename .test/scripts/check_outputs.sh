#!/usr/bin/env bash
# Compare checksums of selected pipeline outputs against a stored manifest
# to detect changes in results between code versions.
#
# Usage: .test/scripts/check_outputs.sh <test> [--update]
#   <test>    name of the test, i.e. .test/config/config_<test>.yaml
#   --update  write the current checksums to .test/expected/<test>.md5
#
# Run from the repository root after the pipeline has finished.
# Only plain text result files are compared, as they are written by the
# pipeline. Compressed files, VCF/BCF files and tool statistics are not
# included because they contain headers with dates, paths or command lines.

set -euo pipefail
shopt -s globstar nullglob

test_name="${1:-}"
mode="${2:-check}"

if [ -z "$test_name" ]; then
    echo "Usage: $0 <test> [--update]" >&2
    exit 2
fi

case "$test_name" in
    mitogenomes)
        patterns=(
            "results/historical/mitogenomes_mapping/stats/ratios_of_*.txt"
        )
        ;;
    gerp)
        patterns=(
            "results/gerp/**/*relative_mutational_load*_table.txt"
        )
        ;;
    mlRho_options)
        patterns=(
            "results/*/mlRho/*/*.mlRho.txt"
        )
        ;;
    pca_roh)
        patterns=(
            "results/*/vcf/*/stats/vcf_merged_missing/*.vcf.stats.txt",
            "results/*/ROH/*FROH_min_2Mb_table.txt"
        )
        ;;
    snpeff)
        patterns=(
            "results/*/snpEff/*.snpEff_variant_impact_table.txt"
        )
        ;;
    userprovided_bams)
        patterns=(
            "results/*/mapping/*/stats/bams_user_provided/*.userprovided.repma.Q30.bam.dpstats.txt",
            "results/*/vcf/*/stats/vcf_merged_missing/*.vcf.stats.txt"
        )
        ;;
    *)
        echo "Error: no output patterns defined for test '$test_name'" >&2
        exit 2
        ;;
esac

# Exclude files based on MapDamage2 rescaling
files=()
for pattern in "${patterns[@]}"; do
    for f in $pattern; do
        case "$(basename "$f")" in
            *.rescaled.*) continue ;;
        esac
        files+=("$f")
    done
done

if [ ${#files[@]} -eq 0 ]; then
    echo "Error: no output files found for test '$test_name'" >&2
    exit 1
fi

expected=".test/expected/${test_name}.md5"
generated="results/checksums/${test_name}.md5"
mkdir -p "$(dirname "$generated")"

printf '%s\n' "${files[@]}" | sort -u | xargs -d '\n' md5sum > "$generated"

if [ "$mode" = "--update" ]; then
    mkdir -p "$(dirname "$expected")"
    cp "$generated" "$expected"
    echo "Wrote $(wc -l < "$expected") checksums to $expected"
    exit 0
fi

if [ ! -f "$expected" ]; then
    echo "::warning::No checksum manifest $expected found. Checksums of this run are in $generated."
    exit 0
fi

if diff -u "$expected" "$generated"; then
    echo "All $(wc -l < "$expected") checksums match $expected"
else
    echo "::error::Outputs of test '$test_name' differ from $expected (lines with - are expected, + are from this run)." \
        "If the change is intended, update the manifest with: $0 $test_name --update"
    exit 1
fi
