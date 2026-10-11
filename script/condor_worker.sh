#!/usr/bin/env bash
set -euo pipefail

# 参数由 analysis.dag 的 VARS 提供；每个节点的 ProcId 都可能是 0。
die() { echo "错误：$*" >&2; exit 1; }
(( $# == 4 )) || die "用法：condor_worker.sh analysis|merge shard result.root expected-count"
mode=$1
shard=$2
result=$3
expected=$4
[[ $mode == analysis || $mode == merge ]] || die "未知的作业模式。"
[[ $shard =~ ^[0-9]+$ && $expected =~ ^[0-9]+$ ]] || die "分片编号/预期数量无效。"

scratch=$(pwd -P)
source "$scratch/job_config.sh"
(( shard < SHARD_COUNT )) || die "分片编号超出范围。"
if [[ $mode == analysis ]]; then
    [[ $result == "part_${shard}.root" ]] || die "分片输出名不匹配。"
else
    [[ $result == "$OUTPUT_NAME" && $expected == "$SHARD_COUNT" ]] || die "合并输出/数量不匹配。"
fi

[[ -r $ROOT_SETUP ]] || die "找不到普通 ROOT 环境：$ROOT_SETUP"
# LCG/thisroot.sh 可能读取未定义变量；只在加载环境时关闭 nounset。
set +u
source "$ROOT_SETUP"
set -u
command -v root >/dev/null || die "加载环境后仍未找到 root。"
echo "Worker: $(hostname); ROOT setup: $ROOT_SETUP; mode: $mode; shard: $shard"
if command -v root-config >/dev/null; then
    root-config --version
fi

tar -xzf "$scratch/analysis.tar.gz" -C "$scratch"
mkdir -p "$scratch/work"
output="$scratch/work/$result"
cd "$scratch/analysis/script"

root_quote() {
    local value=$1
    value=${value//\\/\\\\}
    value=${value//\"/\\\"}
    printf '%s' "$value"
}

output_arg=$(root_quote "$output")
if [[ $mode == analysis ]]; then
    if [[ $MACRO == calculate_c22deltapt ]]; then
        # -1 让入口从分析 JSON 的 c22_delta_pt_output.use_pure 读取模式。
        call=$(printf '%s.cpp("../config/input.json","%s",%s,%s,-1,"../config/analysis.json",%s,%s)' \
            "$MACRO" "$output_arg" "$MAX_FILES" "$MAX_CONFIGS" "$shard" "$SHARD_COUNT")
    else
        call=$(printf '%s.cpp("../config/input.json","%s",%s,%s,"../config/analysis.json",%s,%s)' \
            "$MACRO" "$output_arg" "$MAX_FILES" "$MAX_CONFIGS" "$shard" "$SHARD_COUNT")
    fi
else
    [[ $(wc -l < "$scratch/merge_inputs.txt") -eq $SHARD_COUNT ]] || die "合并清单数量不正确。"
    list_arg=$(root_quote "$scratch/merge_inputs.txt")
    call=$(printf 'merge_results.cpp("%s","%s")' "$output_arg" "$list_arg")
fi

# -n 避免用户 rootlogon.C 加载已安装的 O2Physics 同名类。
if ! root -n -l -b -q "$call" 2>&1 | tee "$scratch/work/root.log"; then
    die "ROOT 运行失败。"
fi
if grep -qE '^Skipped missing input file:|^Error in <|^Fatal in <|\*\*\* Break \*\*\*' "$scratch/work/root.log"; then
    die "发现缺失输入或 ROOT I/O 错误，拒绝上传部分结果。"
fi
if [[ $mode == analysis ]]; then
    loaded=$(grep -c '^Added input file: ' "$scratch/work/root.log" || true)
    [[ $loaded -eq $expected ]] || die "预期 $expected 个输入文件，实际加入 $loaded 个。"
    grep -qE '^Processed [0-9]+ AMPT events; wrote ' "$scratch/work/root.log" || die "分析未正常结束。"
else
    grep -q "^Merged $SHARD_COUNT files into " "$scratch/work/root.log" || die "合并未正常结束。"
fi

[[ -s $output ]] || die "输出为空或不存在。"
task_arg=$(root_quote "$TASK_DIR")
# 重新打开已关闭的文件，拒绝 zombie、恢复过的文件或缺少任务目录的输出。
check_call=$(printf '{ TFile f("%s", "READ"); auto d = f.GetDirectory("%s"); gSystem->Exit(f.IsZombie() || f.TestBit(TFile::kRecovered) || !d || d->GetListOfKeys()->GetSize() == 0 ? 1 : 0); }' \
    "$output_arg" "$task_arg")
if ! root -n -l -b -q -e '#include "../external/o2physics/FlowContainerRoot.h"' -e "$check_call"; then
    die "输出 ROOT 文件检查失败。"
fi

# 只有验证后的文件才放到 transfer_output_files 指定的顶层位置。
mv -- "$output" "$scratch/$result"
echo "Verified output: $result; HTCondor will transfer it to EOS."
