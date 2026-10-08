#!/bin/bash

# 检查是否传入目录参数
if [ -z "$1" ]; then
    echo "❌ 用法: ./mac-check-md5.sh <rawData目录路径>"
    echo "示例: ./mac-check-md5.sh ./data.338F_806R/rawData"
    echo "示例: ./mac-check-md5.sh /Users/xxx/Downloads/MJ_CESHI/data_338F_806R/rawData"
    exit 1
fi

RAW_DIR="$1"

# 检查目录是否存在
if [ ! -d "$RAW_DIR" ]; then
    echo "❌ 错误: 目录不存在 -> $RAW_DIR"
    exit 1
fi

# 进入 rawData 目录
cd "$RAW_DIR" || exit 1

# 查找所有 *md5*.txt 文件（不区分大小写）
LIST_FILE=""
for pattern in "md5.txt" "md5sum.txt" "raw_md5.txt" "clear_md5.txt" "*md5*.txt"; do
    FOUND=$(find . -maxdepth 1 -type f -iname "$pattern" | head -1)
    if [ -n "$FOUND" ]; then
        LIST_FILE="$FOUND"
        break
    fi
done

# 如果还是没找到，报错退出
if [ -z "$LIST_FILE" ] || [ ! -f "$LIST_FILE" ]; then
    echo "❌ 错误: 在 $RAW_DIR 找不到任何 *md5*.txt 校验文件"
    echo "   支持的文件名: md5.txt, md5sum.txt, raw_md5.txt, clear_md5.txt, 或任意 *md5*.txt"
    exit 1
fi

# 去掉开头的 ./
LIST_FILE="${LIST_FILE#./}"

# 先统计总行数（用于进度显示）
# 过滤掉空行和注释行（以#开头的行）
TOTAL_LINES=$(grep -vE '^[[:space:]]*$|^#' "$LIST_FILE" | wc -l | tr -d ' ')

echo "📂 校验目录: $RAW_DIR"
echo "📄 校验清单: $LIST_FILE"
echo "📊 总文件数: $TOTAL_LINES"
echo "================================"

FAIL_COUNT=0
TOTAL_COUNT=0
MISSING_COUNT=0
UNKNOWN_FORMAT_COUNT=0
CURRENT=0

while read -r line; do
    # 跳过空行和注释行
    [ -z "$line" ] && continue
    echo "$line" | grep -qE '^[[:space:]]*#' && continue
    
    CURRENT=$((CURRENT + 1))
    
    # 尝试解析不同的格式
    # 格式1: "MD5 文件名" 或 "MD5 *文件名"
    if echo "$line" | grep -qE '^[a-fA-F0-9]{32}\s+'; then
        hash=$(echo "$line" | awk '{print $1}')
        if echo "$line" | grep -q '\*'; then
            file=$(echo "$line" | awk '{print $3}')
        else
            file=$(echo "$line" | awk '{print $2}')
        fi
    
    # 格式2: "MD5 (文件名) = 校验值" (Mac md5 默认格式)
    elif echo "$line" | grep -qE '^MD5\s*\(.*\)\s*=\s*[a-fA-F0-9]{32}'; then
        file=$(echo "$line" | sed -E 's/^MD5\s*\(([^)]*)\)\s*=\s*[a-fA-F0-9]{32}$/\1/')
        hash=$(echo "$line" | sed -E 's/^MD5\s*\([^)]*\)\s*=\s*([a-fA-F0-9]{32})$/\1/')
    
    # 格式3: "文件名: MD5" (冒号分隔)
    elif echo "$line" | grep -qE '^[^:]+:\s*[a-fA-F0-9]{32}'; then
        file=$(echo "$line" | awk -F': ' '{print $1}')
        hash=$(echo "$line" | awk -F': ' '{print $2}')
    
    # 格式4: "文件名  MD5" (文件名和MD5反过来)
    elif echo "$line" | grep -qE '^[^[:space:]]+\s+[a-fA-F0-9]{32}$'; then
        file=$(echo "$line" | awk '{print $1}')
        hash=$(echo "$line" | awk '{print $2}')
    
    else
        echo "⚠️ [$CURRENT/$TOTAL_LINES] 无法识别的行格式: $line"
        UNKNOWN_FORMAT_COUNT=$((UNKNOWN_FORMAT_COUNT + 1))
        continue
    fi
    
    # 去除可能的 \r 字符（Windows换行符残留）
    file=$(echo "$file" | tr -d '\r')
    hash=$(echo "$hash" | tr -d '\r' | tr '[:upper:]' '[:lower:]')
    
    # 如果文件名为空或哈希值为空，跳过
    if [ -z "$file" ] || [ -z "$hash" ]; then
        echo "⚠️ [$CURRENT/$TOTAL_LINES] 解析失败: $line"
        UNKNOWN_FORMAT_COUNT=$((UNKNOWN_FORMAT_COUNT + 1))
        continue
    fi
    
    TOTAL_COUNT=$((TOTAL_COUNT + 1))
    
    if [ -f "$file" ]; then
        computed=$(md5 -q "$file" | tr '[:upper:]' '[:lower:]')
        if [ "$computed" = "$hash" ]; then
            echo "✅ [$CURRENT/$TOTAL_LINES] $file: OK"
        else
            echo "❌ [$CURRENT/$TOTAL_LINES] $file: FAILED"
            echo "   期望: $hash"
            echo "   实际: $computed"
            FAIL_COUNT=$((FAIL_COUNT + 1))
        fi
    else
        # 尝试在当前目录下查找（如果文件名包含路径，提取basename）
        basename_file=$(basename "$file" 2>/dev/null)
        if [ -n "$basename_file" ] && [ -f "$basename_file" ]; then
            computed=$(md5 -q "$basename_file" | tr '[:upper:]' '[:lower:]')
            if [ "$computed" = "$hash" ]; then
                echo "✅ [$CURRENT/$TOTAL_LINES] $basename_file: OK (文件名不匹配，但内容正确)"
            else
                echo "❌ [$CURRENT/$TOTAL_LINES] $file: 文件不存在 (当前目录下也找不到)"
                echo "   期望: $hash"
                MISSING_COUNT=$((MISSING_COUNT + 1))
                FAIL_COUNT=$((FAIL_COUNT + 1))
            fi
        else
            echo "⚠️ [$CURRENT/$TOTAL_LINES] $file: 文件不存在"
            MISSING_COUNT=$((MISSING_COUNT + 1))
            FAIL_COUNT=$((FAIL_COUNT + 1))
        fi
    fi
done < "$LIST_FILE"

echo "================================"
echo "📊 校验结果:"
echo "   总文件数: $TOTAL_COUNT"
echo "   ✅ 通过: $((TOTAL_COUNT - FAIL_COUNT))"
echo "   ❌ 失败: $FAIL_COUNT"
if [ $MISSING_COUNT -gt 0 ]; then
    echo "   ⚠️ 缺失文件: $MISSING_COUNT"
fi
if [ $UNKNOWN_FORMAT_COUNT -gt 0 ]; then
    echo "   ⚠️ 无法识别的行: $UNKNOWN_FORMAT_COUNT (请检查校验文件格式)"
fi

if [ $FAIL_COUNT -eq 0 ] && [ $UNKNOWN_FORMAT_COUNT -eq 0 ]; then
    echo "🎉 全部通过！数据完整。"
elif [ $FAIL_COUNT -eq 0 ] && [ $UNKNOWN_FORMAT_COUNT -gt 0 ]; then
    echo "⚠️ 文件校验通过，但有 $UNKNOWN_FORMAT_COUNT 行格式无法识别，请检查。"
else
    echo "❌ 有 $FAIL_COUNT 个文件校验失败/缺失，请检查并重新下载。"
fi