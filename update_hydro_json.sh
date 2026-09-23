#!/bin/bash
# update_hydro_json.sh
# 查找 hydro_*_*.root, 提取 id, 统计每个 id 的文件总数, 并生成/更新 JSON

set -euo pipefail

# ============ 可配置项 ============
SEARCH_DIR="/home/huinaibing/hydro_result"
BASE_PATH="/home/huinaibing/hydro_result"
OUTPUT_JSON="/home/huinaibing/git_repo/amptDataAnalysis/config/cent_cfg.json"
STEP=200
# =================================

# 检查搜索目录
[ -d "$SEARCH_DIR" ] || { echo "目录不存在: $SEARCH_DIR" >&2; exit 1; }

# 1. 一次遍历全部文件, 直接按 id 计数
# 避免为每个 id 重复扫描整个搜索目录
declare -A id_count

while IFS= read -r base; do
    if [[ "$base" =~ ^hydro_([0-9]+)_[0-9]+\.root$ ]]; then
        id="${BASH_REMATCH[1]}"
        id_count[$id]=$(( ${id_count[$id]:-0} + 1 ))
    fi
done < <(find "$SEARCH_DIR" -maxdepth 1 -type f -name "hydro_*_*.root" -printf '%f\n')

if [ ${#id_count[@]} -eq 0 ]; then
    echo "未找到匹配文件 hydro_*_*.root" >&2
    echo "[]" > "$OUTPUT_JSON"
    exit 0
fi

mapfile -t ids < <(printf '%s\n' "${!id_count[@]}" | sort -n)

for id in "${ids[@]}"; do
    real=${id_count[$id]}
    id_count[$id]=$(( (real + STEP - 1) / STEP * STEP ))
    echo "id=$id  实际=$real  ->  记录=${id_count[$id]}"
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
        printf '        "n_files": %s\n' "${id_count[$id]}"
        printf '    }' "${id_count[$id]}"
    done
    echo ""
    echo "]"
} > "$OUTPUT_JSON"

echo "已写入: $OUTPUT_JSON  (共 ${#id_count[@]} 个 id)"
