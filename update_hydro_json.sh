#!/bin/bash
# update_hydro_json.sh
# 查找 hydro_*_*.root, 提取每个 id 的文件编号范围, 并生成/更新 JSON

set -euo pipefail

# ============ 可配置项 ============
SEARCH_DIR="/home/huinaibing/hydro_result/"
BASE_PATH="/home/huinaibing/hydro_result/"
OUTPUT_JSON="/home/huinaibing/git_repo/amptDataAnalysis/config/cent_cfg.json"
# =================================

# 检查搜索目录
[ -d "$SEARCH_DIR" ] || { echo "目录不存在: $SEARCH_DIR" >&2; exit 1; }

# 1. 一次遍历全部文件, 记录每个 id 的最小和最大文件编号
# 避免为每个 id 重复扫描整个搜索目录
declare -A id_min id_max

while IFS= read -r base; do
    if [[ "$base" =~ ^hydro_([0-9]+)_([0-9]+)\.root$ ]]; then
        id="${BASH_REMATCH[1]}"
        number=$((10#${BASH_REMATCH[2]}))
        if [[ ! -v id_min[$id] ]] || (( number < id_min[$id] )); then
            id_min[$id]=$number
        fi
        if [[ ! -v id_max[$id] ]] || (( number > id_max[$id] )); then
            id_max[$id]=$number
        fi
    fi
done < <(find "$SEARCH_DIR" -maxdepth 1 -type f -name "hydro_*_*.root" -printf '%f\n')

if [[ -z ${!id_max[*]} ]]; then
    echo "未找到匹配文件 hydro_*_*.root" >&2
    echo "[]" > "$OUTPUT_JSON"
    exit 0
fi

mapfile -t ids < <(printf '%s\n' "${!id_max[@]}" | sort -n)

for id in "${ids[@]}"; do
    span=$(( id_max[$id] - id_min[$id] + 1 ))
    echo "id=$id  编号=${id_min[$id]}..${id_max[$id]}  范围文件数=$span  ->  记录=$(( id_max[$id] + 1 ))"
done

{
    echo "["
    first=1
    for id in "${ids[@]}"; do
        if [ $first -eq 0 ]; then
            echo ","
        fi
        first=0
        printf '    {\n'
        printf '        "path": "%shydro_%s_",\n' "$BASE_PATH" "$id"
        printf '        "bin_val": 0,\n'
        printf '        "n_files": %s\n' "$(( id_max[$id] + 1 ))"
        printf '    }'
    done
    echo ""
    echo "]"
} > "$OUTPUT_JSON"

echo "已写入: $OUTPUT_JSON  (共 ${#id_max[@]} 个 id)"
