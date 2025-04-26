#!/usr/bin/env bash

# 启用严格模式
set -euo pipefail
IFS=$'\n\t'

# 配置路径
PROTEIN_DB=/home/xiaobai/work/Pan_112/Pan_112.fa
BLASTP_DB=/home/xiaobai/work/Pan_112/Pan_112
# 颜色
GREEN='\033[0;32m'
RED='\033[0;31m'
RESET='\033[0m'

print_usage() {
    echo "GenefamilyAna: Run GeneFamily search easily."
    echo ""
    echo "Usage:"
    echo "  $0 -i <Your family fasta file> [-o <genefamily folder name>] [-t <threads>] [-m <i|u>] [-e <BLASTp identity threshold>] [yes/no]"
    echo "    -i <FASTA file>         Input FASTA file"
    echo "    -o <output folder>      Output folder"
    echo "    -t <threads>            Number of threads"
    echo "    -m <mode>               GeneFamily Mode (i=intersection, u=union)"
    echo "    -e <identity threshold> BLASTp identity threshold"
    echo "    [yes/no]                Automatically build tree (optional)"
    echo ""
    echo "Examples:"
    echo "  $0 -i AtSRS_10.fa -o SRS -t 4 -m i -e 75 yes"
    echo "  $0 -i AtSRS_10.fa -m u no"
    echo "  $0 -i AtSRS_10.fa"
    echo "Version: V1.3 (2025-04-26)"
    echo "Contact: 1205630141@qq.com || Author: Ruimiao Li"
    exit 1
}

# 若没有输入任何参数，则打印使用说明
if [ $# -eq 0 ]; then
    print_usage
fi

# 默认参数设定
mode="i"                   # 默认模式为交集(intersection)
threads=1                  # 默认线程数为1
blast_identity_threshold=0 # 默认 BLASTp identity threshold

check_args() {
    # 至少需要 2 个参数：-i 和 FASTA 文件；最多 10 个参数（包括所有可选参数）
    if [[ $# -lt 2 ]] || [[ $# -gt 10 ]]; then
        print_usage
    fi
}

# 在处理任何选项之前调用 check_args 以确保参数数量合理
check_args "$@"

# 判断是否使用新参数形式（检测第一个参数是否以 '-' 开头）
if [[ "$1" == -* ]]; then
    while getopts ":i:o:t:m:e:" opt; do
        case ${opt} in
        i)
            fasta_file="$OPTARG"
            ;;
        o)
            out_folder="$OPTARG"
            ;;
        t)
            threads="$OPTARG"
            ;;
        m)
            mode="$OPTARG"
            ;;
        e)
            blast_identity_threshold="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            print_usage
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            print_usage
            ;;
        esac
    done
fi

# 检查 fasta_file 是否存在
if [[ ! -f "$fasta_file" ]]; then
    echo "Error: $fasta_file does not exist." >&2
    exit 1
fi

# 如果未提供 -o，则使用 fasta_file 的 basename（去除扩展名）作为输出文件夹
if [ -z "${out_folder:-}" ]; then
    out_folder=$(basename "$fasta_file")
    out_folder="${out_folder%.*}"
fi

# 检查额外的 'yes' 或 'no' 参数
shift $((OPTIND - 1))
build_tree_flag="${1:-no}" # 默认为 no

# 打印解析后的参数
echo "FASTA file: $fasta_file"
echo "Protein DB file: $PROTEIN_DB"
echo "Output folder: $out_folder"
echo "Threads: $threads"
echo "GeneFamily Mode (i=intersection, u=union): $mode"
echo "BLASTp identity threshold: ${blast_identity_threshold}"
echo "build_tree_flag: ${build_tree_flag}"
##########################################

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
    local threshold="${3:-0}" # 如果未传递 threshold，则使用默认值 0
    local threads=$4

    local out_file="${family}_blastp.out"
    local blastp_list="blastp.final.${family}.list"

    printf "%b[BLASTp Pipeline]%b\n" "$GREEN" "$RESET"
    if ! blastp -query "$fasta" -db "$BLASTP_DB" -outfmt 7 -evalue 1e-5 -num_threads "$threads" -max_target_seqs 1000000 >"$out_file"; then # -seg yes
        printf "%b[ERROR]%b: BLASTp search failed for %s. Check your files.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi
    printf "BLASTp search for family %b%s%b completed successfully.\n" "$GREEN" "$family" "$RESET"

    if ! grep -v "#" "$out_file" | awk -v thr="$threshold" '$3 > thr {print $2}' | sort -u >"$blastp_list"; then
        printf "%b[ERROR]%b: Filtering BLASTp results failed.\n" "$RED" "$RESET" >&2
        return 1
    fi
    local blastpmethod
    blastpmethod=$(grep -v '^\s*$' "$blastp_list" | wc -l)
    printf "blastp method selected %b%s%b genes.\n" "$GREEN" "$blastpmethod" "$RESET"
    printf "Filtered BLASTp results for family %b%s%b completed successfully.\n" "$GREEN" "$family" "$RESET"
}

run_hmm() {
    local family=$1 # 不与blastp步骤的 $family 冲突
    local mode=$2   # 交集或并集的模式（i 或 u）
    local threads=$3

    local hmm_count # HMM 文件的数量

    printf "%b[HMM Pipeline]%b\n" "$GREEN" "$RESET"

    # 判断当前目录是否包含 HMM 文件
    if ! ls *.hmm 1>/dev/null 2>&1; then
        printf "%b[ERROR]%b: No HMM files found in the current directory ?\n" "$RED" "$RESET" >&2
        return 1
    fi

    # 执行 HMM 搜索
    for hmm_file in *.hmm; do
        hmm=$(basename "$hmm_file" .hmm)

        local out_log="${hmm}.1st.out.log"
        local out_file="${hmm}.1st.out"
        local list_file="${hmm}_1st_qua_id.txt"

        # 执行 HMM 搜索
        if ! hmmsearch --cpu $threads --noali --cut_tc -o "$out_log" --domtblout "$out_file" "$hmm_file" "$PROTEIN_DB"; then
            printf "%b[ERROR]%b: HMM search failed for %s.\n" "$RED" "$RESET" "$hmm_file" >&2
            return 1
        fi

        # 提取高置信度的序列
        if ! grep -v "#" "$out_file" | awk '($7 + 0) < 1E-5' | cut -f1 -d " " | sort -u >"$list_file"; then
            printf "%b[ERROR]%b: Extracting confidence list failed for %s.\n" "$RED" "$RESET" "$hmm_file" >&2
            return 1
        fi
        hmmsearch_count=$(grep -v '^\s*$' "$list_file" | wc -l)
        printf "HMM search for %b%s%b completed successfully, selected ${hmmsearch_count} members.\n" "$GREEN" "$hmm_file" "$RESET"
    done

    # 计算 HMM 文件数量
    hmm_count=$(find . -maxdepth 1 -type f -name "*_1st_qua_id.txt" | wc -l)

    # 根据传入的模式选择交集或并集
    if [[ "$mode" == "i" || "$mode" == "I" ]]; then
        # 处理交集
        if [[ $hmm_count -eq 1 ]]; then
            printf "Only %b%s%b HMM file found. No need to intersect, copying directly.\n" "$GREEN" "$hmm_count" "$RESET"
            cp ./*_1st_qua_id.txt "${family}_hmm_1st.list"
        elif [[ $hmm_count -gt 1 ]]; then
            printf "%b%s%b HMM files found. Extracting intersection...\n" "$GREEN" "$hmm_count" "$RESET"
            if ! awk -v hmm_count="$hmm_count" '{ a[$1]++ } END { for (i in a) if (a[i] == hmm_count) print i }' ./*_1st_qua_id.txt | sort >"${family}_hmm_1st.list"; then
                printf "%b[ERROR]%b: Failed to extract the HMM intersection.\n" "$RED" "$RESET" >&2
                return 1
            fi
            printf "Successfully extracted the intersection of %b%s%b HMM files.\n" "$GREEN" "$hmm_count" "$RESET"
        fi
    elif [[ "$mode" == "u" || "$mode" == "U" ]]; then
        # 处理并集
        printf "Extracting union of %b%s%b HMM files...\n" "$GREEN" "$hmm_count" "$RESET"
        if ! cat ./*_1st_qua_id.txt | sort | uniq >"${family}_hmm_1st.list"; then
            printf "%b[ERROR]%b: Failed to extract union.\n" "$RED" "$RESET" >&2
            return 1
        fi
        printf "Successfully extracted the union of %b%s%b HMM files.\n" "$GREEN" "$hmm_count" "$RESET"
    else
        # 如果选择无效，默认选择交集
        echo "Invalid input. Defaulting to intersection."
        if [[ $hmm_count -eq 1 ]]; then
            cp ./*_1st_qua_id.txt "${family}_hmm_1st.list"
        elif [[ $hmm_count -gt 1 ]]; then
            if ! awk -v hmm_count="$hmm_count" '{ a[$1]++ } END { for (i in a) if (a[i] == hmm_count) print i }' ./*_1st_qua_id.txt | sort >"${family}_hmm_1st.list"; then
                echo "%b[ERROR]%b: Failed to extract intersection." >&2
                return 1
            fi
        fi
        printf "Intersection mode has been applied by default.\n"
    fi
}

extract_sequences_4_rebuild_HMM() {
    local family=$1
    local threads=$2

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
    if ! seqkit grep -f "$source_file" "$PROTEIN_DB" -o "$output_file" -j $threads; then
        printf "%s[ERROR]%s: Extracting sequences failed for '%s'.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi

    printf "Extracted sequences for family %b%s%b successfully.\n" "$GREEN" "$family" "$RESET"

}

multi_align_seq() {
    local family=$1
    local threads=$2

    local source_file="${family}_hmm_1st.fa"
    local output_file="${family}_hmm_1st.fa.aln"

    if ! mafft --auto --quiet --thread $threads "$source_file" >"$output_file"; then
        printf "%b[ERROR]%b: mafft multiple sequence alignment failed." "$RED" "$RESET" >&2
        return 1
    fi
    echo "mafft multiple sequence alignment for HMM rebuilt completed ."
}

rebuild_hmm() {
    local family=$1
    local threads=$2

    local rebuild_hmm="${family}_rebuild.hmm"
    local rebuild_log="${family}_rebuild.hmm.log"
    local out_log="${family}.rebuild.out.log"
    local out_file="${family}.rebuild.out"
    local list_file="HMM.final.${family}.list"
    #
    if ! hmmbuild --cpu $threads --amino -o "$rebuild_log" "$rebuild_hmm" "${family}_hmm_1st.fa.aln"; then
        printf "%b[ERROR]%b: rebulit HMM failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi
    printf "rebuilt HMM completed.\n"
    #
    if ! hmmsearch --cpu $threads --noali --domtblout "$out_file" -o "$out_log" "$rebuild_hmm" "$PROTEIN_DB"; then
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
    hmmmethod=$(grep -v '^\s*$' "$list_file" | wc -l) # 不算空行
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
    # 创建目标目录，如果不存在
    if [[ ! -d "$destination_dir" ]]; then
        mkdir -p "$destination_dir" || {
            printf "Error: Failed to create directory %s\n" "$destination_dir" >&2
            return 1
        }
    fi

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
    local target_dir="./${family}/RawData"

    # 创建目标目录，如果不存在
    if ! mkdir -p "$target_dir"; then
        printf "%b[ERROR]%b: Failed to create directory %s.\n" "$RED" "$RESET" "$target_dir" >&2
        return 1
    fi

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

    # 移动文件
    if ! mv "$fasta" "$target_dir/" || ! mv $hmm_file "$target_dir/"; then
        printf "%b[ERROR]%b: Moving RawData files failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi

    printf "%b[SUCCESS]%b: RawData for family %s moved successfully to %s.\n" "$GREEN" "$RESET" "$family" "$target_dir"
}

# 树构建函数
build_tree() {
    local family=$1
    local aln_file="./${family}/Tree/common.${family}.fa.aln"
    local tree_file="./${family}/Tree/${family}.nwk"

    mkdir -p "./${family}/Tree"

    # 对齐序列
    if ! mafft --auto --quiet --thread 8 "./${family}/common.${family}.fa" >"$aln_file"; then
        printf "%b[ERROR]%b: Alignment failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi

    # 构建树
    if ! fasttree -log "$tree_file".log "$aln_file" >"$tree_file"; then
        printf "%b[ERROR]%b: Tree construction failed for %s.\n" "$RED" "$RESET" "$family" >&2
        return 1
    fi
    printf "Tree data for family %b%s%b created successfully.\n" "$GREEN" "$family" "$RESET"
}
# 此步基础在基因家族成员鉴定完毕之后，因此将在目的文件夹中进行函数操作。

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
    done <./110SampleID.txt

    # 计算众数
    local mode
    mode=$(printf "%s\n" "${counts[@]}" | sort -n | uniq -c | sort -nr | awk 'NR==1 {print $2}')

    printf "Sample\tType\tGenes\n" >"$result_file"

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

        printf "%s\t%s\t%s\n" "$sample_id" "$sample_type" "$gene_list" >>"$result_file"
    done
}

main() {
    # 获取传递给 main 的参数
    local fasta="$1"
    local family="$2"
    local threads="$3"
    local mode="$4"
    local threshold="$5"
    local build_tree_flag="$6" # 新增参数，用来接收 yes/no 选项

    # 预先检查参数
    check_enviro

    # 获取脚本名称
    SCRIPT_NAME=$(basename "$0")
    printf "[%s] All the necessary commands are verified, running...\n" "$SCRIPT_NAME"

    # 创建目标目录
    mkdir -p "$family"

    # 运行 BLAST 比对
    run_blastp "$family" "$fasta" "$threshold" "$threads"

    # 运行 HMM
    run_hmm "$family" "$mode" "$threads"

    # 其他分析步骤
    extract_sequences_4_rebuild_HMM "$family" "$threads"
    multi_align_seq "$family" "$threads"
    rebuild_hmm "$family" "$threads"

    # 交集方法
    common_method "$family"

    # 移动文件
    move_files "$family" "./$family"
    move_rawData "$family" "$fasta"

    # 树构建
    # Normalize build_tree_flag input (ignore case and extra spaces)
    build_tree_flag_normalized=$(echo "$build_tree_flag" | tr '[:upper:]' '[:lower:]' | xargs)

    # Check if build_tree_flag is "yes" or "no" after normalization
    if [[ "$build_tree_flag_normalized" =~ ^(yes|y)$ ]]; then
        echo "Automatically building tree for family $out_folder..."
        build_tree "$out_folder"
    elif [[ "$build_tree_flag_normalized" =~ ^(no|n)$ ]]; then
        echo "Skipping tree building for family $out_folder."
    else
        echo "%b[ERROR]%b: Invalid build_tree_flag value: '$build_tree_flag'. Please use 'yes' or 'no'." >&2
        return 1
    fi

    # 统计泛基因组成员
    count_Pan_members "$family"
}

main "$fasta_file" "$out_folder" "$threads" "$mode" "$blast_identity_threshold" "$build_tree_flag"