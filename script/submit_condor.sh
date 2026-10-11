#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)

# ============ 可配置项：修改这里后直接运行，无需命令行参数 ============
MACRO="calculate_v2ptrho"  # 与 run_parallel.sh 支持相同的七个入口
INPUT_JSON="$REPO_DIR/config/cent_cfg_ampt_eos.json"
ANALYSIS_JSON="$REPO_DIR/config/config.json"
OUTPUT_NAME="myAnalysisResultV2PtRho.root"
EOS_INPUT_ENDPOINT="root://eosuser.cern.ch"  # 转换 JSON 中的 /eos/ 前缀
EOS_OUTPUT_DIR="root://eosuser.cern.ch//eos/user/q/qiuyu/ampt/condor"
ROOT_SETUP="/cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh"
SHARD_COUNT=1000
MAX_JOBS=1000                 # DAG 中同时排队/运行的节点上限
MAX_FILES=-1               # 每个原始 JSON entry 的文件上限；-1 表示全部
MAX_CONFIGS=-1             # 原始 JSON entry 数量上限；-1 表示全部
RETRIES=2                  # 每个节点失败后最多重试两次
SUBMIT_DIR="$SCRIPT_DIR/condor_runs"  # CERN 上应位于 AFS，不要放在 EOS
SUBMIT=true                # false：只生成文件；true：生成后提交 DAG
# ==================================================================

die() { echo "错误：$*" >&2; exit 1; }
(( $# == 0 )) || die "请修改脚本开头的变量，不需要传入参数。"
command -v python3 >/dev/null || die "未找到 python3。"
[[ $SHARD_COUNT =~ ^[1-9][0-9]*$ && $MAX_JOBS =~ ^[1-9][0-9]*$ ]] || die "SHARD_COUNT/MAX_JOBS 必须是正整数。"
[[ $MAX_FILES =~ ^(-1|0|[1-9][0-9]*)$ && $MAX_CONFIGS =~ ^(-1|0|[1-9][0-9]*)$ ]] || die "数量上限必须为 -1 或非负整数。"
[[ $RETRIES =~ ^[0-9]+$ ]] || die "RETRIES 必须为非负整数。"
[[ $SUBMIT == true || $SUBMIT == false ]] || die "SUBMIT 必须为 true 或 false。"
[[ -f $INPUT_JSON && -f $ANALYSIS_JSON ]] || die "找不到输入/分析 JSON；请先生成 EOS 输入清单。"

if [[ $SUBMIT == true ]]; then
    command -v condor_submit_dag >/dev/null || die "未找到 condor_submit_dag，请在 CERN 提交节点运行。"
    command -v klist >/dev/null && klist -s || die "没有有效 Kerberos ticket；请先运行 kinit。"
fi

mkdir -p -- "$SUBMIT_DIR"
run_dir=$(mktemp -d "$SUBMIT_DIR/run_$(date +%Y%m%d_%H%M%S)_XXXXXX")

# 打包源码和两份配置，不传输原始 ROOT 数据，也不在提交节点跑分析。
python3 - "$REPO_DIR" "$run_dir" "$MACRO" "$INPUT_JSON" "$ANALYSIS_JSON" \
    "$OUTPUT_NAME" "$EOS_INPUT_ENDPOINT" "$EOS_OUTPUT_DIR" "$ROOT_SETUP" \
    "$SHARD_COUNT" "$MAX_FILES" "$MAX_CONFIGS" "$RETRIES" <<'PY'
import io
import json
import math
from pathlib import Path
import re
import shlex
import shutil
import sys
import tarfile

(repo, run, macro, input_json, analysis_json, output_name, endpoint,
 eos_dir, root_setup, shards, max_files, max_configs, retries) = sys.argv[1:]
repo, run = Path(repo).resolve(), Path(run).resolve()
shards, max_files, max_configs = map(int, (shards, max_files, max_configs))
task_dirs = {
    "calculate_v2ptrho": "pid-flow-pt-corr",
    "calculate_c22deltapt": "pid-flow-pt-corr",
    "calculate_pidptcorrelations": "pid-flow-pt-corr",
    "calculate_flowGfwOmegaXi": "flow-gfw-omega-xi",
    "calculate_flowPbpbPikp": "flow-pbpb-pikp_TPC_TOF",
    "calculate_ptspectra": "pt-spectra",
    "calculate_pid_v2_rp_gap_variance": "pid-v2-rp-gap-moments",
}

def check(condition, message):
    if not condition:
        sys.exit("错误：" + message)

check(macro in task_dirs, "不支持的 MACRO。")
check(not str(run).startswith("/eos/"), "提交文件和日志必须放在 AFS/本地目录。")
# DAG and submit macros use unquoted simple paths; reject ambiguous characters.
check(re.fullmatch(r"/[A-Za-z0-9_./-]+", str(run)), "SUBMIT_DIR 不能包含空格或特殊字符。")
check(re.fullmatch(r"[A-Za-z0-9_.-]+\.root", output_name), "OUTPUT_NAME 必须是简单的 .root 文件名。")
check(not re.fullmatch(r"part_[0-9]+\.root", output_name), "最终文件名不能与分片文件名相同。")
check(re.fullmatch(r"root://[A-Za-z0-9_.:-]+", endpoint), "EOS_INPUT_ENDPOINT 格式无效。")
check(re.fullmatch(r"root://[A-Za-z0-9_.:-]+//eos/[A-Za-z0-9_./-]+", eos_dir), "EOS_OUTPUT_DIR 必须是 EOS XRootD 目录 URI。")
check(root_setup.startswith("/") and "\n" not in root_setup, "ROOT_SETUP 必须是 worker 可访问的绝对路径。")
check(max(shards, max_files, max_configs) <= 2147483647, "数量参数超过 ROOT int 范围。")

with open(input_json) as source:
    inputs = json.load(source)
check(isinstance(inputs, list), "输入 JSON 顶层必须是数组。")
for entry in inputs:
    check(isinstance(entry, dict), "输入 JSON entry 必须是对象。")
    prefix, count, bin_val = (entry.get(key) for key in ("path", "n_files", "bin_val"))
    check(isinstance(prefix, str) and prefix and "\n" not in prefix and "\r" not in prefix, "path 必须是非空单行字符串。")
    check(type(count) is int and 0 <= count <= 2147483647, "n_files 必须是非负 int。")
    check(type(bin_val) in (int, float) and math.isfinite(bin_val), "bin_val 必须是有限数值。")
    if prefix.startswith("/eos/"):
        entry["path"] = endpoint + "/" + prefix
    else:
        check(prefix.startswith(("root://", "/")), "path 必须为绝对路径或 XRootD URI。")

# 保留所有 entry 的原始顺序，bootstrap 使用原始 config index/file number。
selected = inputs if max_configs < 0 else inputs[:max_configs]
counts = [entry["n_files"] if max_files < 0 else min(entry["n_files"], max_files)
          for entry in selected]
expected = [sum(len(range(shard, count, shards)) for count in counts)
            for shard in range(shards)]
check(sum(expected) > 0, "所选输入没有文件。")
with open(analysis_json) as source:
    analysis_bytes = source.read().encode()
json.loads(analysis_bytes)  # 打包前检查 JSON 语法。
input_bytes = (json.dumps(inputs, indent=2, ensure_ascii=False) + "\n").encode()

with tarfile.open(run / "analysis.tar.gz", "w:gz") as archive:
    for header in sorted(repo.glob("*.h")):
        archive.add(header, arcname="analysis/" + header.name)
    for folder in ("dataFrame", "external"):
        archive.add(repo / folder, arcname="analysis/" + folder)
    for source in sorted((repo / "script").glob("*.cpp")):
        archive.add(source, arcname="analysis/script/" + source.name)
    for name, data in (("input.json", input_bytes), ("analysis.json", analysis_bytes)):
        info = tarfile.TarInfo("analysis/config/" + name)
        info.size, info.mode = len(data), 0o644
        archive.addfile(info, io.BytesIO(data))

settings = {
    "MACRO": macro, "TASK_DIR": task_dirs[macro], "OUTPUT_NAME": output_name,
    "ROOT_SETUP": root_setup, "SHARD_COUNT": str(shards),
    "MAX_FILES": str(max_files), "MAX_CONFIGS": str(max_configs),
}
(run / "job_config.sh").write_text("".join(key + "=" + shlex.quote(value) + "\n"
                                           for key, value in settings.items()))
eos_output = eos_dir.rstrip("/") + "/" + run.name + "/"
(run / "merge_inputs.txt").write_text("".join(eos_output + f"part_{shard}.root\n"
                                             for shard in range(shards)))
for name in ("run.sub", "condor_worker.sh"):
    shutil.copy2(repo / "script" / name, run / name)
(run / "condor_worker.sh").chmod(0o755)
(run / "logs").mkdir()

lines = []
parents = []
for shard in range(shards):
    node = f"SHARD_{shard}"
    parents.append(node)
    lines += [f"JOB {node} run.sub",
              f'VARS {node} mode="analysis" shard="{shard}" result="part_{shard}.root" expected="{expected[shard]}" eos_output="{eos_output}"',
              f"RETRY {node} {retries}"]
lines += ["JOB MERGE run.sub",
          f'VARS MERGE mode="merge" shard="0" result="{output_name}" expected="{shards}" eos_output="{eos_output}"',
          f"RETRY MERGE {retries}",
          "PARENT " + " ".join(parents) + " CHILD MERGE"]
(run / "analysis.dag").write_text("\n".join(lines) + "\n")
print(f"作业文件：{run}")
print(f"输入文件总数：{sum(expected)}；分片数：{shards}；空分片数：{expected.count(0)}")
print(f"最终输出：{eos_output}{output_name}")
PY

cd -- "$run_dir"
if [[ $SUBMIT == true ]]; then
    condor_submit_dag -maxjobs "$MAX_JOBS" analysis.dag
else
    printf '仅生成文件。提交命令：cd %q && condor_submit_dag -maxjobs %s analysis.dag\n' "$PWD" "$MAX_JOBS"
fi
