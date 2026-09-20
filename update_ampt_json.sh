#!/bin/bash
# update_ampt_json.sh
# 查找 ampt_*_0.root,提取 id,统计每个 id 的文件总数,并生成/更新 JSON

set -euo pipefail

# ============ 可配置项 ============
SEARCH_DIR="/home/huinaibing/new_ampt_version/"                                              # 搜索目录
BASE_PATH="/home/huinaibing/new_ampt_version/"              # JSON 中 path 前缀
OUTPUT_JSON="/home/huinaibing/git_repo/amptDataAnalysis/config/cent_cfg.json"                             # 输出 JSON 路径
STEP=200                                                    # 向上取整的步长(比如 100
# =================================

# 检查搜索目录
[ -d "$SEARCH_DIR" ] || { echo "目录不存在: $SEARCH_DIR" >&2; exit 1; }

round_up() {
    local n=$1
    echo $(( (n + STEP - 1) / STEP * STEP ))
}


# 1. find . -name "ampt_*_0.root" 拿到种子文件,提取 id
# 使用关联数组去重,顺便记录每个 id 的文件总数
declare -A id_count

while IFS= read -r file; do
    base=$(basename "$file")
    id=$(echo "$base" | sed -E 's/^ampt_([0-9]+)_0\.root$/\1/')
    [ -z "$id" ] && continue

    if [ -z "${id_count[$id]:-}" ]; then
        real=$(find "$SEARCH_DIR" -maxdepth 1 -type f -name "ampt_${id}_*.root" | wc -l)
        id_count[$id]=$(round_up "$real")
        echo "id=$id  实际=$real  ->  记录=${id_count[$id]}"
    fi
done < <(find "$SEARCH_DIR" -maxdepth 1 -type f -name "ampt_*_0.root" | sort)

if [ ${#id_count[@]} -eq 0 ]; then
    echo "未找到匹配文件 ampt_*_0.root" >&2
    echo "[]" > "$OUTPUT_JSON"
    exit 0
fi

{
    echo "["
    first=1
    for id in $(printf '%s\n' "${!id_count[@]}" | sort -n); do
        if [ $first -eq 0 ]; then
            echo ","
        fi
        first=0
        printf '    {\n'
        printf '        "path": "%sampt_%s_",\n' "$BASE_PATH" "$id"
        printf '        "bin_val": 0,\n'
        printf '        "n_files": %s\n' "${id_count[$id]}"
        printf '    }'
    done
    echo ""
    echo "]"
} > "$OUTPUT_JSON"

echo "已写入: $OUTPUT_JSON  (共 ${#id_count[@]} 个 id)"
