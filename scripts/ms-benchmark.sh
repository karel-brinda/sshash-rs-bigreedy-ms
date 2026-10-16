#!/usr/bin/env bash
# Compare SSHash-rs indexes built from two representations of one k-mer set:
# a plain SPSS (unitigs / simplitigs / eulertigs, as given) versus a BiGreedy
# masked superstring of the same k-mers (computed here with KmerCamel).
#
# Reports, for each representation: the index size with its component
# breakdown (strings / offsets / MPHF / skew index), point-lookup timings
# (positive and negative), and, when a query file is given, streaming-query
# timings. The masked-superstring index is also verified against its input
# (`check --masked-superstring`: k-mer ids must be the ranks of the
# uppercase positions, which fails on any non-min-one mask).
#
# Usage:
#   scripts/ms-benchmark.sh <spss.fa[.gz]> <k> <m> [query.fa/fq[.gz]]
#
# Requirements: kmercamel in PATH (https://github.com/ondrejsladky/kmercamel);
# the sshash binary is taken from $SSHASH, or target/release/sshash.

set -euo pipefail

if [[ $# -lt 3 ]]; then
    grep '^#' "$0" | sed 's/^# \{0,1\}//' | head -16
    exit 1
fi

INPUT=$1
K=$2
M=$3
QUERY=${4:-}

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SSHASH=${SSHASH:-"$SCRIPT_DIR/../target/release/sshash"}
if [[ ! -x "$SSHASH" ]]; then
    echo "error: sshash binary not found at $SSHASH (build with: cargo build --release)" >&2
    exit 1
fi
command -v kmercamel >/dev/null || { echo "error: kmercamel not in PATH" >&2; exit 1; }

WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/sshash-ms-bench.XXXXXX")
trap 'rm -rf "$WORKDIR"' EXIT
MSFA="$WORKDIR/ms.msfa"

echo "## Computing the BiGreedy masked superstring (kmercamel, k=$K)"
kmercamel compute -k "$K" -o "$MSFA" "$INPUT"
echo

breakdown() {
    perl -pe 's/\e\[[0-9;]*m//g' \
        | grep "sshash_lib::dictionary" \
        | grep -E "total index size|\[bits/kmer\]" \
        | sed -E 's/^.*sshash_lib::dictionary:[[:space:]]*//'
}

echo "## Building from the plain SPSS: $INPUT"
"$SSHASH" build --input "$INPUT" -k "$K" -m "$M" --output "$WORKDIR/spss" 2>&1 | breakdown
echo
echo "## Building from the masked superstring: $MSFA"
"$SSHASH" build --input "$MSFA" -k "$K" -m "$M" --masked-superstring --output "$WORKDIR/ms" 2>&1 | breakdown
echo

echo "## Verifying the masked-superstring index (min-one id bijection)"
"$SSHASH" check --index "$WORKDIR/ms" --input "$MSFA" --masked-superstring 2>&1 | grep -E "Total k-mers|Failures|CHECK"
echo

for rep in spss ms; do
    echo "## Point lookups ($rep)"
    "$SSHASH" bench --index "$WORKDIR/$rep" 2>&1 | grep -E "positive lookup|negative lookup"
    echo
done

if [[ -n "$QUERY" ]]; then
    for rep in spss ms; do
        echo "## Streaming queries ($rep) over $QUERY"
        "$SSHASH" stream-bench --index "$WORKDIR/$rep" --query "$QUERY" 2>&1 \
            | grep -E "num_kmers|found_kmers|time_per_kmer|extension_ratio"
        echo
    done
fi
