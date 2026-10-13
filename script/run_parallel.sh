#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <macro> <workers> [input-json] [output-root] [max-files] [max-configs] [analysis-json]" >&2
    echo "Macros: calculate_v2ptrho, calculate_c22deltapt, calculate_pidptcorrelations, calculate_flowGfwOmegaXi, calculate_flowPbpbPikp, calculate_ptspectra, calculate_pid_v2_rp_gap_variance, calculate_pid_charged_sc" >&2
    exit 2
}

(( $# >= 2 && $# <= 7 )) || usage
macro=$1
workers=$2
default_input=../config/cent_cfg.json
case "$macro" in
    calculate_v2ptrho) default_output=myAnalysisResultV2PtRho.root ;;
    calculate_c22deltapt) default_output=myAnalysisResultC22DeltaPt_pt04.root ;;
    calculate_pidptcorrelations) default_output=myAnalysisResultPidPtCorrelations.root ;;
    calculate_flowGfwOmegaXi) default_output=myAnalysisResultFlowGfwOmegaXi.root ;;
    calculate_flowPbpbPikp) default_output=myAnalysisResultFlowPbpbPikp.root ;;
    calculate_ptspectra) default_output=myAnalysisResultPtSpectra.root ;;
    calculate_pid_v2_rp_gap_variance) default_output=pid_v2_rp_gap_moments.root; default_input=../config/cent_cfg_ampt.json ;;
    calculate_pid_charged_sc) default_output=pid_charged_sc_moments.root; default_input=../config/cent_cfg_ampt.json ;;
    *) usage ;;
esac
[[ $workers =~ ^[1-9][0-9]*$ ]] || usage

input=${3:-$default_input}
output=${4:-$default_output}
max_files=${5:--1}
max_configs=${6:--1}
analysis=${7:-../config/config.json}
[[ $max_files =~ ^-?[0-9]+$ && $max_configs =~ ^-?[0-9]+$ ]] || usage
command -v root >/dev/null || { echo "ROOT is not available; load ROOT first" >&2; exit 1; }

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
        # -1 让入口从分析 JSON 的 c22_delta_pt_output.use_pure 读取模式。
        call=$(printf '%s.cpp("%s","%s",%s,%s,-1,"%s",%d,%d)' \
            "$macro" "$input" "$part_arg" "$max_files" "$max_configs" "$analysis" "$shard" "$workers")
    else
        call=$(printf '%s.cpp("%s","%s",%s,%s,"%s",%d,%d)' \
            "$macro" "$input" "$part_arg" "$max_files" "$max_configs" "$analysis" "$shard" "$workers")
    fi
    root -n -l -b -q "$call" >"$part_dir/part_${shard}.log" 2>&1 &
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

printf '%s\n' "${parts[@]}" > "$part_dir/inputs.txt"
output_arg=$(root_quote "$output")
input_list_arg=$(root_quote "$part_dir/inputs.txt")
merge_call=$(printf 'merge_results.cpp("%s","%s")' "$output_arg" "$input_list_arg")
if ! root -n -l -b -q "$merge_call"; then
    echo "Merge failed; logs and partial files are in $part_dir" >&2
    exit 1
fi
if [[ $macro == calculate_pid_charged_sc ]]; then
    # Reject stale/partial results that cannot supply the charged baseline.
    check_call=$(printf '{ TFile f("%s", "READ"); auto d=f.GetDirectory("pid-charged-sc/Charged"); gSystem->Exit(f.IsZombie() || !d || !d->Get("pTwoCharged") || !d->Get("pFourCharged") || !d->Get("pTwoChargedByGroup") || !d->Get("pFourChargedByGroup") ? 1 : 0); }' "$output_arg")
    if ! root -n -l -b -q -e "$check_call"; then
        echo "Missing charged SC baseline components; logs and partial files are in $part_dir" >&2
        exit 1
    fi
fi
echo "Merged output: $output"
rm -rf "$part_dir"
