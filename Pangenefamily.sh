#!/bin/bash

# set -euo pipefail # 捕获潜在的错误

# 配置路径和颜色
PROTEIN_DB=/home/xiaobai/work/Pan110_xiaomi/Pan110_xiaomi.fa # ~/work/PanSi_110/110pro.primarytranscript.fa
BLASTP_DB=/home/xiaobai/work/Pan110_xiaomi/Pan110_xiaomi
# color
GREEN='\033[0;32m'
RED='\033[0;31m'
RESET='\033[0m'

# %s 用于输出普通字符串。
# %b 用于处理字符串中的转义字符并输出相应的效果。

print_usage() {
    printf "GenefamilyAna: Run GeneFamily search easily.\n"
    printf "Usage: bash %s <genefamily folder name> <Your family fasta file> [BLASTp identity threshold]\n" "$0"
    printf "E.g.: bash %s SRS AtSRS_10.fa\n" "$0"
    printf "Version: V1.3 (2024-11-19)\n"
    printf "Contact me: 1205630141@qq.com || Author: Ruimiao Li\n"
    exit 1
}

check_args() {
    if [[ $# -lt 2 ]] || [[ $# -gt 3 ]]; then
        print_usage
    fi
}

check_enviro() {
    local blastp hmmscan mafft seqkit fasttree

    blastp=$(command -v blastp)
    if [[ -z "$blastp" || ! -x "$blastp" ]]; then
        echo "Could not find blastp in the server"
        return 1
    fi

    hmmscan=$(command -v hmmscan)
    if [[ -z "$hmmscan" || ! -x "$hmmscan" ]]; then
        echo "Could not find hmmscan in the server"
        return 1
    fi

    mafft=$(command -v mafft)
    if [[ -z "$mafft" || ! -x "$mafft" ]]; then
        echo "Could not find mafft in the server"
        return 1
    fi

    seqkit=$(command -v seqkit)
    if [[ -z "$seqkit" || ! -x "$seqkit" ]]; then
        echo "Could not find seqkit in the server"
        return 1
    fi

    fasttree=$(command -v fasttree)
    if [[ -z "$fasttree" || ! -x "$fasttree" ]]; then
        echo "Could not find fasttree in the server"
        return 1
    fi
}

run_blastp() {
    local family=$1
    local fasta=$2

    local out_file="${family}_blastp.out"
    local blastp_list="blastp.final.${family}.list"

    printf "%b[1]%b: BLASTp Pipeline\n" "$GREEN" "$RESET"
    if ! blastp -query "$fasta" -db "$BLASTP_DB" -outfmt 7 -evalue 1e-5 -num_threads 8 -seg yes -max_target_seqs 1000000 >"$out_file"; then
        printf "%b[ERROR]%b: BLASTp search failed for %s. Check your files.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi
    printf "BLASTp search for family %b%s%b completed successfully.\n" "$GREEN" "$family" "$RESET"

    local threshold="${3:-0}" # variable="${param:-default}" 如果 param 未设置或为空，则使用 default 作为值；如果 param 有值，则直接使用该值。
    if ! grep -v "#" "$out_file" | awk -v thr="$threshold" '$3 > thr {print $2}' | sort -u >"$blastp_list"; then
        printf "%b[ERROR]%b: Filtering BLASTp results failed.\n" "$RED" "$RESET" >&2
        return 1
    fi
    local blastpmethod
    blastpmethod=$(wc -l <"$blastp_list")
    printf "blastp method selected %b%s%b genes.\n" "$GREEN" "$blastpmethod" "$RESET"
    printf "Filtered BLASTp results for family %b%s%b completed successfully.\n" "$GREEN" "$family" "$RESET"
}

run_hmm() {
    local family=$1 # 不与blastp步骤的 $family 冲突

    printf "%b[2]%b: HMM Pipeline\n" "$GREEN" "$RESET"

    # 循环处理每个HMM文件
    for hmm_file in *.hmm; do
        hmm=$(basename "$hmm_file" .hmm)

        local out_log="${hmm}.1st.out.log"
        local out_file="${hmm}.1st.out"
        local list_file="${hmm}_1st_qua_id.txt"

        # 执行 HMM 搜索
        if ! hmmsearch --cpu 8 --noali --cut_tc -o "$out_log" --domtblout "$out_file" "$hmm_file" "$PROTEIN_DB"; then
            printf "%b[ERROR]%b: HMM search failed for %s.\n" "$RED" "$RESET" "$hmm_file" >&2
            return 1
        fi

        # 提取具有高置信度的序列
        if ! grep -v "#" "$out_file" | awk '($7 + 0) < 1E-5' | cut -f1 -d " " | sort -u >"$list_file"; then
            printf "%b[ERROR]%b: Extracting confidence list failed for %s.\n" "$RED" "$RESET" "$hmm_file" >&2
            return 1
        fi

        # 成功处理该文件的消息
        printf "HMM search for %b%s%b completed successfully.\n" "$GREEN" "$hmm_file" "$RESET"
    done

    # 计算符合条件的HMM文件数量
    local hmm_count

    hmm_count=$(find . -maxdepth 1 -type f -name "*_1st_qua_id.txt" | wc -l)

    # 根据HMM文件数量处理交集提取
    if [[ $hmm_count -eq 1 ]]; then
        printf "Only %b%s%b HMM file found, so no intersection needed.\n" "$GREEN" "$hmm_count" "$RESET"
        cp ./*_1st_qua_id.txt "${family}_hmm_1st.list" # mv
    elif [[ $hmm_count -gt 1 ]]; then
        printf " %b%s%b HMM files found, getting intersection ...\n" "$GREEN" "$hmm_count" "$RESET"
        # 提取多个HMM文件交集
        if ! awk -v hmm_count="$hmm_count" '{a[$1]++} END{for (i in a) if (a[i] == hmm_count) print i}' ./*_1st_qua_id.txt | sort >"${family}_hmm_1st.list"; then # awk -v hmm_count="$hmm_count" 'FNR==1{a[$1]++} END{for (i in a) if(a[i]==hmm_count) print i}' ./*_1st_qua_id.txt | sort >"${family}_hmm_1st.list"
            printf "%b[ERROR]%b: Intersection of %s HMM confidence lists failed.\n" "$RED" "$RESET" "$hmm_count" >&2
            return 1
        fi
        printf "Successfully extracted the intersection of the confidence sequence lists from %b%s%b HMM files for %b%s%b.\n" "$GREEN" "$hmm_count" "$RESET" "$GREEN" "$family" "$RESET"
    fi
    ###

}

extract_sequences_4_rebuild_HMM() {
    local family=$1
    local source_file="${family}_hmm_1st.list"
    local output_file="${family}_hmm_1st.fa"
    # if [[ ! -f "$PROTEIN_DB" ]]; then
    # printf "%s[ERROR]%s: Protein database file '%s' not found.\n" "$RED" "$RESET" "$PROTEIN_DB" >&2
    # return 1  # 或 exit 1 如果在脚本主流程
    # local faidx="$PROTEIN_DB.fai"
    # if [ ! -x "$faidx" ]; then
    #     seqkit faidx "$PROTEIN_DB"
    # fi

    # 然后进行grep操作
    if ! seqkit grep -f "$source_file" "$PROTEIN_DB" -o "$output_file" -j 4; then
        printf "%s[ERROR]%s: Extracting sequences failed for '%s'.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi

    printf "Extracted sequences for family %b%s%b successfully.\n" "$GREEN" "$family" "$RESET"

}

multi_align_seq() {
    local family=$1
    local source_file="${family}_hmm_1st.fa"
    local output_file="${family}_hmm_1st.fa.aln"

    if ! mafft --auto --quiet --thread 4 "$source_file" >"$output_file"; then
        printf "%b[ERROR]%b: mafft multiple sequence alignment failed." "$RED" "$RESET" >&2
        return 1
    fi
    echo "mafft multiple sequence alignment for HMM rebuilt completed ."
}

rebuild_hmm() {
    local family=$1
    local rebuild_hmm="${family}_rebuild.hmm"
    local rebuild_log="${family}_rebuild.hmm.log"
    local out_log="${family}.rebuild.out.log"
    local out_file="${family}.rebuild.out"
    local list_file="HMM.final.${family}.list"
    #
    if ! hmmbuild --cpu 8 --amino -o "$rebuild_log" "$rebuild_hmm" "${family}_hmm_1st.fa.aln"; then
        printf "%b[ERROR]%b: rebulit HMM failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi
    printf "rebuilt HMM completed.\n"
    #
    if ! hmmsearch --cpu 8 --noali --domtblout "$out_file" -o "$out_log" "$rebuild_hmm" "$PROTEIN_DB"; then
        printf "%b[ERROR]%b: Extracting sequences failed for %s.\n" "$RED" "$RESET" "$family" >&2
        rm "$rebuild_hmm"
        return 1
    fi
    # 提取具有高置信度的序列
    if ! grep -v "#" "$out_file" | awk '($7 + 0) < 1E-5' | cut -f1 -d " " | sort -u >"$list_file"; then
        printf "%b[ERROR]%b: Extracting confidence list failed for %s.\n" "$RED" "$RESET" "$hmm_file" >&2
        return 1
    fi
    # 统计选定的HMM方法数量
    local hmmmethod
    hmmmethod=$(wc -l <"$list_file")
    printf "HMM method selected %b%d%b genes.\n" "$GREEN" "$hmmmethod" "$RESET"
}

common_method() {
    local family=$1

    local common_list="common.${family}.list"
    local genefamily_fasta="common.${family}.fa"

    comm -12 <(sort "blastp.final.${family}.list") <(sort "HMM.final.${family}.list") >"$common_list"
    common_method=$(wc -l <"$common_list")
    printf "Common sequence candidates for %b%s%b include %b%s%b genes.\n" "$GREEN" "$family" "$RESET" "$GREEN" "$common_method" "$RESET"
    seqkit grep -f "$common_list" $PROTEIN_DB -o "$genefamily_fasta"
}

move_files() {
    local family="$1"
    local destination_dir="$2"
    local pattern_matched

    # 定义文件清单
    local file_patterns=(
        "${family}_blastp.out"
        "blastp.final.${family}.list"

        "*.1st.out.log"
        "*.1st.out"
        "*_1st_qua_id.txt"

        "${family}_hmm_1st.list"
        "${family}_hmm_1st.fa"
        "${family}_hmm_1st.fa.aln"
        "${family}_rebuild.hmm.log"
        "${family}_rebuild.hmm"
        "${family}.rebuild.out.log"
        "${family}.rebuild.out"
        "HMM.final.${family}.list"

        "common.${family}.list"
        "common.${family}.fa"
    )

    mkdir -p "$destination_dir"
    # 处理文件匹配
    for pattern in "${file_patterns[@]}"; do
        # 检查是否有匹配的文件
        if pattern_matched=$(find . -maxdepth 1 -type f -name "$pattern" 2>/dev/null); then
            if [[ -n "$pattern_matched" ]]; then
                # 移动匹配到的所有文件
                mv $pattern_matched "$destination_dir" || {
                    printf "Error: Failed to move files matching %s to %s\n" "$pattern" "$destination_dir" >&2
                    return 1
                }
            else
                printf "Warning: No files matching pattern %s\n" "$pattern" >&2
            fi
        else
            printf "Warning: No files matching pattern %s\n" "$pattern" >&2
        fi
    done
    printf "Files moved to %s successfully.\n" "$destination_dir"

}

move_rawData() {
    local family="$1"
    local fasta="$2"
    local hmm_file
    local target_dir="./${family}/RawData"

    # 查找 HMM 文件并验证存在性
    hmm_file=$(find . -maxdepth 1 -type f -name "*.hmm" 2>/dev/null)
    if [[ -z "$hmm_file" ]]; then
        printf "%b[ERROR]%b: No HMM files found for family %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi

    # 验证 FASTA 文件存在性
    if [[ ! -e "$fasta" ]]; then
        printf "%b[ERROR]%b: Fasta file %s does not exist.\n" "$RED" "$RESET" "$fasta" >&2
        return 1
    fi

    # 创建目标目录
    if ! mkdir -p "$target_dir"; then
        printf "%b[ERROR]%b: Failed to create directory %s.\n" "$RED" "$RESET" "$target_dir" >&2
        return 1
    fi

    # 移动文件
    if ! mv "$fasta" "$target_dir/" || ! mv $hmm_file "$target_dir/"; then
        printf "%b[ERROR]%b: Moving RawData files failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi

    printf "%b[SUCCESS]%b: RawData for family %s moved successfully to %s.\n" "$GREEN" "$RESET" "$family" "$target_dir"
    return 0
}

# 此步基础在基因家族成员鉴定完毕之后，因此将在目的文件夹中进行函数操作。
build_tree() {
    local family=$1

    local aln_file="./${family}/Tree/common.${family}.fa.aln"
    local tree_file="./${family}/Tree/${family}.nwk"

    mkdir -p "./${family}/Tree"
    if ! mafft --auto --quiet --thread 8 "./${family}/common.${family}.fa" >"$aln_file"; then
        printf "%b[ERROR]%b: Alignment failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi

    if ! fasttree -log "$tree_file".log "$aln_file" >"$tree_file"; then
        printf "%b[ERROR]%b: Tree construction failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi
    printf "Tree data for family %b%s%b created successfully.\n" "$GREEN" "$family" "$RESET"
}

ask_to_build_tree() {
    local choice
    local family=$1

    while true; do
        printf "Do you want to build the tree for family %s? (y[Y]/n[N]): " "$family"
        read -r choice

        if [[ "$choice" =~ ^[Yy]$ ]]; then
            if ! build_tree "$family"; then
                printf "Error: Failed to build the tree for family %s.\n" "$family" >&2
                return 1
            fi
            break
        elif [[ "$choice" =~ ^[Nn]$ ]]; then
            printf "Skipping tree building for family %s.\n" "$family"
            break
        elif [[ -z "$choice" ]]; then
            printf "No input detected. Skipping tree building for family %s.\n" "$family"
            break
        else
            printf "Invalid input. Please enter 'y' for yes or 'n' for no.\n"
        fi
    done
}

count_Pan_members() {
    local family=$1
    local common_list="./${family}/common.${family}.list" 
    local result_file="./${family}/${family}_Pan_result.xls"

    # 初始化用于存储样本 ID 和计数的数组
    local sample_ids=()
    local counts=()

    # 检查依赖文件是否存在
    if [[ ! -f ./110SampleID.txt ]]; then
        printf "Error: 110SampleID.txt not found.\n" >&2
        return 1
    fi

    if [[ ! -f "$common_list" ]]; then
        printf "Error: %s not found.\n" "$common_list" >&2
        return 1
    fi

    # 读取样本 ID 并统计计数
    while IFS= read -r sample_id; do
        if [[ -z "$sample_id" ]]; then
            continue
        fi

        local count
        if ! count=$(grep -c "${sample_id}_" "$common_list" 2>/dev/null); then
            count=0
        fi

        sample_ids+=("$sample_id")
        counts+=("$count")
    done < ./110SampleID.txt

    # 计算众数
    local mode
    mode=$(printf "%s\n" "${counts[@]}" | sort -n | uniq -c | sort -nr | awk 'NR==1 {print $2}')

    printf "Sample\tType\tGenes\n" > "$result_file"

    # 遍历样本 ID 和计数，分类输出结果
    for ((i = 0; i < ${#sample_ids[@]}; i++)); do
        local sample_id="${sample_ids[$i]}"
        local count="${counts[$i]}"
        local sample_type=""
        local gene_list

        # 判断分类
        if [[ $count -gt $mode ]]; then
            sample_type="more than $mode"
        elif [[ $count -lt $mode ]]; then
            sample_type="less than $mode"
        else
            sample_type="$mode members"
        fi

        # 提取基因列表
        if ! gene_list=$(grep "${sample_id}_" "$common_list" | tr '\n' '\t' | sed 's/\t$//'); then
            gene_list="None"
        fi

        printf "%s\t%s\t%s\n" "$sample_id" "$sample_type" "$gene_list" >> "$result_file"
    done
}

main() {
    # 预先检查参数
    check_args "$@"
    check_enviro

    local family=$1
    local fasta=$2
    local threshold="${3:-0}"

    # 创建目标目录
    mkdir -p "$family"
    
    # 运行 BLAST 比对
    run_blastp "$family" "$fasta" "$threshold"

    # 运行 HMM
    run_hmm "$family"
    extract_sequences_4_rebuild_HMM "$family"
    multi_align_seq "$family"
    rebuild_hmm "$family"

    # 交集方法
    common_method "$family"

    # 移动文件
    move_files "$family" "./$family"
    move_rawData "$family" "$fasta"

    # 可选的树构建
    ask_to_build_tree "$family"

    # 统计泛基因组成员
    count_Pan_members "$family"
}

# 执行 main 函数
main "$@"