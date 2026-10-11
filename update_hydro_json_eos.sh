#!/usr/bin/env bash
# update_hydro_json_eos.sh
# 在 CERN EOS 服务器上扫描 hydro particles 目录并生成 hydro 输入 JSON。
#
# 该脚本故意不使用 find。EOS 目录中的文件数量很大时，挂载路径上的
# find 可能报 "Argument list too long"；EOS 客户端绕过挂载层枚举目录。

set -euo pipefail

# ============ 可配置项 ============
export EOS_MGM_URL="root://eosuser.cern.ch"              # CERN 用户 EOS 服务，避免连接 localhost
SEARCH_DIR="/eos/user/q/qiuyu/hydro/result/PbPb5360/particles/"  # 搜索目录
BASE_PATH="/eos/user/q/qiuyu/hydro/result/PbPb5360/particles/"   # JSON 中 path 前缀
OUTPUT_JSON="./config/cent_cfg_hydro_eos.json"           # 输出 JSON 路径；从项目根目录运行
# =================================

BASE_PATH="${BASE_PATH%/}/"

command -v eos >/dev/null 2>&1 || {
    echo "错误：服务器上未找到 eos 命令，请先加载 EOS/AliEn 环境。" >&2
    exit 1
}

mkdir -p "$(dirname "$OUTPUT_JSON")"
declare -A id_min id_max

# 先检查列目录是否成功；进程替换无法把 eos 的退出状态传给 while。
# -F 在目录名末尾添加 /，避免把形似 ROOT 文件的目录计入结果。
listing=$(mktemp)
trap 'rm -f -- "$listing"' EXIT
if ! eos ls -F "$SEARCH_DIR" > "$listing"; then
    echo "EOS 列目录失败，未更新 JSON: $SEARCH_DIR" >&2
    exit 1
fi

while IFS= read -r base; do
    [[ "$base" =~ ^hydro_([0-9]+)_([0-9]+)\.root$ ]] || continue
    id="${BASH_REMATCH[1]}"
    number=$((10#${BASH_REMATCH[2]}))
    if [[ ! -v id_min[$id] ]] || (( number < id_min[$id] )); then
        id_min[$id]=$number
    fi
    if [[ ! -v id_max[$id] ]] || (( number > id_max[$id] )); then
        id_max[$id]=$number
    fi
done < "$listing"

if [[ -z ${!id_max[*]} ]]; then
    echo "未找到匹配文件 hydro_*_*.root: $SEARCH_DIR" >&2
    exit 1
fi

mapfile -t ids < <(printf '%s\n' "${!id_max[@]}" | sort -n)
{
    echo '['
    first=1
    for id in "${ids[@]}"; do
        ((first)) || echo ','
        first=0
        printf '    {\n'
        printf '        "path": "%shydro_%s_",\n' "$BASE_PATH" "$id"
        printf '        "bin_val": 0,\n'
        printf '        "n_files": %s\n' "$((id_max[$id] + 1))"
        printf '    }'
    done
    echo
    echo ']'
} > "$OUTPUT_JSON"

echo "已写入: $OUTPUT_JSON  (共 ${#id_max[@]} 个 id)"
