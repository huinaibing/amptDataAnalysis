#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <macro> <workers> [input-json] [output-root] [max-files] [max-configs] [analysis-json] [use-pure-for-c22]" >&2
    exit 2
}

(( $# >= 2 && $# <= 8 )) || usage
macro=$1
workers=$2
case "$macro" in
    calculate_v2ptrho) default_output=myAnalysisResultV2PtRho.root ;;
    calculate_c22deltapt) default_output=myAnalysisResultC22DeltaPt_pt04.root ;;
    calculate_pidptcorrelations) default_output=myAnalysisResultPidPtCorrelations.root ;;
    calculate_flowGfwOmegaXi) default_output=myAnalysisResultFlowGfwOmegaXi.root ;;
    *) usage ;;
esac
[[ $workers =~ ^[1-9][0-9]*$ ]] || usage

input=${3:-../config/cent_cfg.json}
output=${4:-$default_output}
max_files=${5:--1}
max_configs=${6:--1}
analysis=${7:-../config/config.json}
use_pure=${8:--1}
[[ $max_files =~ ^-?[0-9]+$ && $max_configs =~ ^-?[0-9]+$ && $use_pure =~ ^-?[0-9]+$ ]] || usage
[[ $macro == calculate_c22deltapt || $# -lt 8 ]] || usage
command -v root >/dev/null || { echo "ROOT is not available; load O2Physics first" >&2; exit 1; }
command -v hadd >/dev/null || { echo "hadd is not available; load O2Physics first" >&2; exit 1; }

cd "$(dirname "$0")"
part_dir=$(mktemp -d "${TMPDIR:-/tmp}/ampt-shards.XXXXXX")
parts=()
pids=()

root_quote() {
    local value=$1
    value=${value//\\/\\\\}
    value=${value//\"/\\\"}
    printf '%s' "$value"
}

input=$(root_quote "$input")
analysis=$(root_quote "$analysis")
for (( shard = 0; shard < workers; ++shard )); do
    part="$part_dir/part_${shard}.root"
    part_arg=$(root_quote "$part")
    parts+=("$part")
    if [[ $macro == calculate_c22deltapt ]]; then
        call=$(printf '%s.cpp("%s","%s",%s,%s,%s,"%s",%d,%d)' \
            "$macro" "$input" "$part_arg" "$max_files" "$max_configs" "$use_pure" "$analysis" "$shard" "$workers")
    else
        call=$(printf '%s.cpp("%s","%s",%s,%s,"%s",%d,%d)' \
            "$macro" "$input" "$part_arg" "$max_files" "$max_configs" "$analysis" "$shard" "$workers")
    fi
    root -l -b -q "$call" >"$part_dir/part_${shard}.log" 2>&1 &
    pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        failed=1
    fi
done
for part in "${parts[@]}"; do
    if [[ ! -s $part ]]; then
        failed=1
    fi
done
if (( failed )); then
    echo "A shard failed; logs and partial files are in $part_dir" >&2
    exit 1
fi

hadd -f "$output" "${parts[@]}"
echo "Merged output: $output"
echo "Shard logs and files: $part_dir"
