#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# 汉坦病毒单样本预处理 + 组装流程
# 流程: fastp → PRINSEQ++ → bowtie2(去rRNA) → MEGAHIT → Prodigal
#
# 用法:
#   bash scripts/qc.sh -s <SAMPLE> [-c config.yaml] [--dry-run]
#   bash scripts/qc.sh -c config.yaml [--dry-run]          (自动发现 data_dir 下所有样本)
#
# 输出:
#   默认写入 <project_dir>/results，可在 config.yaml 的 paths.results_dir 中指定
#   <results_dir>/json/<SAMPLE>.qc.json   — 结构化结果（供下游 skill 衔接）
#   <results_dir>/log/<SAMPLE>/           — 每步独立日志
#   屏幕无输出，静默运行
#
# 可选清理:
#   config.yaml 设置 cleanup.intermediates: true 时，MEGAHIT 完成后删除
#   前三步中间 fastq（fastp 的 html/json 质控报告保留）
# ============================================================

# ---------- 路径初始化 ----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG="$(cd "${SCRIPT_DIR}/.." && pwd)/config.yaml"
PYTHON_BIN="${PYTHON:-python3}"

SAMPLE=""
CONFIG_YAML="${CONFIG_YAML:-${DEFAULT_CONFIG}}"
DRY_RUN="false"
COMPLETED_STEPS=""

# ---------- 帮助信息 ----------
usage() {
    cat <<'EOF'
用法:
  bash scripts/qc.sh -s <样本名> [-c 配置文件] [--dry-run]
  bash scripts/qc.sh <样本名>
  bash scripts/qc.sh -c 配置文件 [--dry-run]          (自动发现 data_dir 下所有样本)
EOF
}

# ---------- JSON 错误输出（写入文件，不打印到屏幕）----------
plain_json_error() {
    local sample="${1:-unknown}" message="$2" code="${3:-1}"
    local err_json="${JSON_DIR:-/tmp}/${sample}.qc.json"
    mkdir -p "$(dirname "${err_json}")" 2>/dev/null || true
    printf '{"%s":{"sample":"%s","status":"failed","message":"%s","exit_code":%s}}\n' \
        "${sample}" "${sample}" "${message//\"/\\\"}" "${code}" > "${err_json}"
    exit "${code}"
}

# ---------- 参数解析 ----------
while [ "$#" -gt 0 ]; do
    case "$1" in
        -s|--sample) SAMPLE="${2:-}"; shift 2 ;;
        -c|--config) CONFIG_YAML="${2:-}"; shift 2 ;;
        --dry-run)   DRY_RUN="true"; shift ;;
        -h|--help)   usage; exit 0 ;;
        -*)
            usage >&2
            plain_json_error "${SAMPLE:-unknown}" "未知选项: $1" 2
            ;;
        *)
            if [ -z "${SAMPLE}" ]; then SAMPLE="$1"; shift
            else usage >&2; plain_json_error "${SAMPLE}" "多余参数: $1" 2
            fi
            ;;
    esac
done

# ---------- 前置检查 ----------
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || plain_json_error "${SAMPLE:-unknown}" "需要 python3 来解析 config.yaml" 127
[ -f "${CONFIG_YAML}" ] || plain_json_error "${SAMPLE:-unknown}" "找不到配置文件: ${CONFIG_YAML}" 2

# ============================================================
# 配置加载：从 config.yaml 读取所有参数，导出为 shell 变量
# ============================================================
load_config() {
    "${PYTHON_BIN}" - "$CONFIG_YAML" "$SAMPLE" <<'PY'
import os, shlex, sys

cfg_path, sample = sys.argv[1], sys.argv[2]

def scalar(value):
    value = value.strip()
    if not value:
        return ""
    if (value[0:1], value[-1:]) in [("'", "'"), ('"', '"')]:
        return value[1:-1]
    return value

def read_simple_yaml(path):
    cfg, section = {}, None
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip()
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            indent = len(line) - len(line.lstrip(" "))
            text = line.strip()
            if ":" not in text:
                continue
            key, value = text.split(":", 1)
            key, value = key.strip(), value.split(" #", 1)[0].strip()
            if indent == 0:
                if value == "":
                    cfg[key] = {}
                    section = key
                else:
                    cfg[key] = scalar(value)
                    section = None
            elif section:
                cfg.setdefault(section, {})[key] = scalar(value)
    return cfg

def get(section, key, default=""):
    data = cfg.get(section, {})
    return data.get(key, default) if isinstance(data, dict) else default

def emit(name, value):
    print(f"{name}={shlex.quote(str(value))}")

cfg = read_simple_yaml(cfg_path)
project_dir = cfg.get("project_dir")
if not project_dir:
    raise SystemExit("ERROR: config.yaml 中必须设置 project_dir")

paths = cfg.get("paths", {}) if isinstance(cfg.get("paths"), dict) else {}
results_dir = paths.get("results_dir") or os.path.join(project_dir, "results")

# --- 目录路径 ---
emit("PROJECT_DIR", project_dir)
emit("DATA_DIR", paths.get("data_dir") or os.path.join(project_dir, "data"))
emit("RESULTS_DIR", results_dir)
emit("FASTP_DIR", paths.get("fastp_dir") or os.path.join(results_dir, "01_fastp"))
emit("PRINSEQ_DIR", paths.get("prinseq_dir") or os.path.join(results_dir, "02_prinseq"))
emit("CLEAN_DIR", paths.get("clean_dir") or os.path.join(results_dir, "03_clean_non_rRNA"))
emit("MEGAHIT_DIR", paths.get("megahit_dir") or os.path.join(results_dir, "04_megahit"))
emit("PRODIGAL_DIR", paths.get("prodigal_dir") or os.path.join(results_dir, "05_prodigal"))
emit("LOG_DIR", paths.get("log_dir") or os.path.join(results_dir, "log", sample))
emit("JSON_DIR", paths.get("json_dir") or os.path.join(results_dir, "json"))

# --- 输入发现 ---
emit("R1_SUFFIXES", get("input", "r1_suffixes", "_r1.fq.gz,_R1.fq.gz,_1.fastq.gz"))
emit("R2_SUFFIXES", get("input", "r2_suffixes", "_r2.fq.gz,_R2.fq.gz,_2.fastq.gz"))

# --- 线程 ---
emit("THREADS_CFG", cfg.get("threads", "auto"))

# --- fastp ---
fp = lambda k, d: get("fastp", k, d)
emit("FASTP_DISABLE_ADAPTER", fp("disable_adapter_trimming", "true"))
emit("FASTP_DETECT_ADAPTER", fp("detect_adapter_for_pe", "false"))
emit("FASTP_CUT_FRONT", fp("cut_front", "true"))
emit("FASTP_CUT_TAIL", fp("cut_tail", "true"))
emit("FASTP_CUT_RIGHT", fp("cut_right", "true"))
w = fp("cut_window_size", "5")
emit("FASTP_FRONT_WINDOW", fp("cut_front_window_size", w))
emit("FASTP_TAIL_WINDOW", fp("cut_tail_window_size", w))
emit("FASTP_RIGHT_WINDOW", fp("cut_right_window_size", w))
q = fp("cut_mean_quality", "20")
emit("FASTP_FRONT_QUAL", fp("cut_front_mean_quality", q))
emit("FASTP_TAIL_QUAL", fp("cut_tail_mean_quality", q))
emit("FASTP_RIGHT_QUAL", fp("cut_right_mean_quality", q))
emit("FASTP_QUALIFIED_PHRED", fp("qualified_quality_phred", "20"))
emit("FASTP_UNQUALIFIED_PERCENT", fp("unqualified_percent_limit", "20"))
emit("FASTP_LENGTH", fp("length_required", "50"))
emit("FASTP_N_LIMIT", fp("n_base_limit", "0"))
emit("FASTP_EXTRA_ARGS", fp("extra_args", ""))

# --- PRINSEQ++ ---
ps = lambda k, d: get("prinseq", k, d)
emit("PRINSEQ_ENTROPY", ps("lc_entropy", "0.5"))
emit("PRINSEQ_DUST", ps("lc_dust", "0.5"))
emit("PRINSEQ_EXTRA_ARGS", ps("extra_args", ""))

# --- bowtie2（替换 SortMeRNA）---
bt2 = lambda k, d: get("bowtie2", k, d)
emit("BOWTIE2_REF_DIR", bt2("ref_dir", os.path.join(project_dir, "DB", "sortmerna")))
emit("BOWTIE2_INDEX", os.path.join(bt2("ref_dir", os.path.join(project_dir, "DB", "sortmerna")), bt2("index_prefix", "rrna_index")))
emit("BOWTIE2_PRESET", bt2("preset", "very-sensitive"))
emit("BOWTIE2_EXTRA_ARGS", bt2("extra_args", ""))

# --- MEGAHIT ---
mh = lambda k, d: get("megahit", k, d)
emit("MEGAHIT_MIN_CONTIG_LEN", mh("min_contig_len", "200"))
emit("MEGAHIT_OVERWRITE", mh("overwrite", "true"))
emit("MEGAHIT_EXTRA_ARGS", mh("extra_args", ""))

# --- Prodigal ---
pd = lambda k, d: get("prodigal", k, d)
emit("PRODIGAL_MODE", pd("mode", "meta"))
emit("PRODIGAL_MIN_AA", pd("min_protein_aa", "200"))
emit("PRODIGAL_EXTRA_ARGS", pd("extra_args", ""))

# --- 中间文件清理（MEGAHIT 完成后删除前三步 fastq；默认关闭）---
cl = lambda k, d: get("cleanup", k, d)
emit("CLEANUP_INTERMEDIATES", cl("intermediates", "false"))
PY
}

CONFIG_ENV="$(load_config)" || plain_json_error "${SAMPLE:-unknown}" "解析 config.yaml 失败: ${CONFIG_YAML}" 2
eval "${CONFIG_ENV}"

# ============================================================
# 线程与并行度计算（提前到这里，供多样本并行调度使用）
# ============================================================

# 自动检测可用线程数（读 config.yaml 的 threads，auto 则自动探测）
detect_threads() {
    if [ "${THREADS_CFG}" != "auto" ] && [ -n "${THREADS_CFG}" ]; then
        printf '%s\n' "${THREADS_CFG}"
    elif command -v nproc >/dev/null 2>&1; then
        nproc
    elif command -v sysctl >/dev/null 2>&1; then
        sysctl -n hw.ncpu
    else
        printf '4\n'
    fi
}

THREADS="$(detect_threads)"

# 每个进程（软件）的线程数：fastp 单进程上限 16，其它软件与之对齐，
# 保证多样本并行时总核数不超过配置值
PER_PROC_THREADS=16
# 并行样本数 = 总线程数 / 每进程线程数（向下取整，至少 1）
# 预防总线程数不是 16 的倍数：余数核闲置，绝不超订
if [ "${THREADS}" -gt "${PER_PROC_THREADS}" ]; then
    PARALLEL=$(( THREADS / PER_PROC_THREADS ))
else
    PARALLEL=1
    PER_PROC_THREADS="${THREADS}"
fi
[ "${PARALLEL}" -ge 1 ] || PARALLEL=1

# ============================================================
# 自动发现样本（未指定 -s 时扫描 data_dir）
# ============================================================
if [ -z "${SAMPLE}" ]; then
    SAMPLES_FOUND=""
    IFS=',' read -r -a _R1_ARR <<< "${R1_SUFFIXES}"
    for _sfx in "${_R1_ARR[@]}"; do
        for _f in "${DATA_DIR}"/*"${_sfx}"; do
            [ -f "${_f}" ] || continue
            _bn="$(basename "${_f}" "${_sfx}")"
            SAMPLES_FOUND="${SAMPLES_FOUND}${SAMPLES_FOUND:+ }${_bn}"
        done
    done
    # 去重
    SAMPLES_FOUND="$(printf '%s\n' ${SAMPLES_FOUND} | sort -u | tr '\n' ' ')"
    SAMPLES_FOUND="${SAMPLES_FOUND% }"

    if [ -z "${SAMPLES_FOUND}" ]; then
        usage >&2
        plain_json_error "unknown" "未指定 -s 且在 ${DATA_DIR} 中未自动发现任何样本（后缀: ${R1_SUFFIXES}）" 2
    fi

    printf '[%s] 自动发现 %s 个样本，并行 %s 路（每样本每软件 %s 线程）\n' \
        "$(date '+%F %T')" "${SAMPLES_FOUND}" "${PARALLEL}" "${PER_PROC_THREADS}" >&2
    EXIT_TOTAL=0
    RUNNING=0
    for _sp in ${SAMPLES_FOUND}; do
        if [ "${DRY_RUN}" = "true" ]; then
            bash "${BASH_SOURCE[0]}" -s "${_sp}" -c "${CONFIG_YAML}" --dry-run &
        else
            bash "${BASH_SOURCE[0]}" -s "${_sp}" -c "${CONFIG_YAML}" &
        fi
        RUNNING=$(( RUNNING + 1 ))
        if [ "${RUNNING}" -ge "${PARALLEL}" ]; then
            if wait -n; then :; else EXIT_TOTAL=1; fi
            RUNNING=$(( RUNNING - 1 ))
        fi
    done
    while [ "${RUNNING}" -gt 0 ]; do
        if wait -n; then :; else EXIT_TOTAL=1; fi
        RUNNING=$(( RUNNING - 1 ))
    done
    exit ${EXIT_TOTAL}
fi

# ============================================================
# 工具函数
# ============================================================

RESULT_JSON="${JSON_DIR}/${SAMPLE}.qc.json"

# 一次性创建所有输出目录
mkdir -p "${FASTP_DIR}" "${PRINSEQ_DIR}" "${CLEAN_DIR}" \
         "${MEGAHIT_DIR}" "${PRODIGAL_DIR}" "${LOG_DIR}" "${JSON_DIR}"

# ---------- JSON 结果写入（只写文件，不打印到屏幕）----------
write_json() {
    local status="$1" failed_step="${2:-}" message="${3:-}" exit_code="${4:-0}"
    STATUS="${status}" FAILED_STEP="${failed_step}" MESSAGE="${message}" EXIT_CODE="${exit_code}" \
    SAMPLE="${SAMPLE}" THREADS="${THREADS}" COMPLETED_STEPS="${COMPLETED_STEPS}" \
    LOG_DIR="${LOG_DIR}" RESULT_JSON="${RESULT_JSON}" \
    FASTP_R1="${FASTP_R1:-}" FASTP_R2="${FASTP_R2:-}" \
    PRINSEQ_R1="${PRINSEQ_R1:-}" PRINSEQ_R2="${PRINSEQ_R2:-}" \
    CLEAN_R1="${CLEAN_R1:-}" CLEAN_R2="${CLEAN_R2:-}" \
    CONTIGS_FASTA="${CONTIGS_FASTA:-}" \
    PROTEIN_FAA="${PROTEIN_FAA:-}" FILTERED_FAA="${FILTERED_FAA:-}" \
    "${PYTHON_BIN}" - <<'PY'
import json, os

sample = os.environ["SAMPLE"]
completed = [x for x in os.environ.get("COMPLETED_STEPS", "").split(",") if x]
failed = os.environ.get("FAILED_STEP", "")
steps = {name: {"status": "ok", "log": os.path.join(os.environ["LOG_DIR"], f"{name}.log")}
         for name in completed}
if failed:
    steps[failed] = {"status": "failed", "log": os.path.join(os.environ["LOG_DIR"], f"{failed}.log")}

outputs = {
    "fastp_r1": os.environ.get("FASTP_R1", ""),
    "fastp_r2": os.environ.get("FASTP_R2", ""),

    "prinseq_r1": os.environ.get("PRINSEQ_R1", ""),
    "prinseq_r2": os.environ.get("PRINSEQ_R2", ""),
    "clean_non_rrna_r1": os.environ.get("CLEAN_R1", ""),
    "clean_non_rrna_r2": os.environ.get("CLEAN_R2", ""),
    "contigs_fasta": os.environ.get("CONTIGS_FASTA", ""),
    "prodigal_faa": os.environ.get("PROTEIN_FAA", ""),
    "prodigal_min_aa_faa": os.environ.get("FILTERED_FAA", ""),
}

payload = {
    sample: {
        "sample": sample,
        "status": os.environ["STATUS"],
        "failed_step": failed,
        "message": os.environ.get("MESSAGE", ""),
        "exit_code": int(os.environ.get("EXIT_CODE", "0")),
        "threads": int(os.environ.get("THREADS", "0") or 0),
        "log_dir": os.environ["LOG_DIR"],
        "result_json": os.environ["RESULT_JSON"],
        "steps": steps,
        "outputs": {k: v for k, v in outputs.items() if v},
    }
}

text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
with open(os.environ["RESULT_JSON"], "w", encoding="utf-8") as f:
    f.write(text + "\n")
# 不打印到 stdout，保持静默
PY
}

# 失败时调用：写入错误 JSON 并退出
fail_step() {
    local step="$1" code="$2" message="$3"
    write_json "failed" "${step}" "${message}" "${code}"
    exit "${code}"
}

# 写入管道日志
log_msg() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" >> "${LOG_DIR}/pipeline.log"
}

# 检查必需工具是否存在
require_tool() {
    local tool="$1"
    command -v "${tool}" >/dev/null 2>&1 || fail_step "preflight" 127 "未找到必需工具: ${tool}"
}

# 追加额外参数（空格分割）
append_extra_args() {
    local extra="$1"
    [ -n "${extra}" ] || return 0
    read -r -a EXTRA_ARGV <<< "${extra}"
    CMD+=("${EXTRA_ARGV[@]}")
}

# 执行一个步骤：记录日志 → 运行命令 → 检查退出码
run_step() {
    local step="$1"; shift
    local log_file="${LOG_DIR}/${step}.log"
    log_msg "START ${step}"

    # 记录命令到日志
    {
        printf '[%s] CMD:' "$(date '+%F %T')"
        printf ' %q' "$@"
        printf '\n'
    } > "${log_file}"

    if [ "${DRY_RUN}" = "true" ]; then
        printf '[%s] DRY_RUN 跳过执行\n' "$(date '+%F %T')" >> "${log_file}"
    else
        set +e
        "$@" >> "${log_file}" 2>&1
        local code=$?
        set -e
        if [ "${code}" -ne 0 ]; then
            log_msg "FAIL ${step} exit=${code}; 日志: ${log_file}"
            fail_step "${step}" "${code}" "${step} 执行失败; 详见 ${log_file}"
        fi
    fi

    COMPLETED_STEPS="${COMPLETED_STEPS:+${COMPLETED_STEPS},}${step}"
    log_msg "DONE ${step}"
}

# 在数据目录中查找样本的双端 FASTQ 文件
find_reads() {
    local r1_suffix r2_suffix candidate_r1 candidate_r2
    IFS=',' read -r -a R1_ARR <<< "${R1_SUFFIXES}"
    IFS=',' read -r -a R2_ARR <<< "${R2_SUFFIXES}"
    [ "${#R1_ARR[@]}" -eq "${#R2_ARR[@]}" ] || \
        fail_step "input" 2 "r1_suffixes 和 r2_suffixes 条目数必须一致"

    for idx in "${!R1_ARR[@]}"; do
        r1_suffix="${R1_ARR[$idx]}"
        r2_suffix="${R2_ARR[$idx]}"
        candidate_r1="${DATA_DIR}/${SAMPLE}${r1_suffix}"
        candidate_r2="${DATA_DIR}/${SAMPLE}${r2_suffix}"
        if [ -f "${candidate_r1}" ] && [ -f "${candidate_r2}" ]; then
            RAW_R1="${candidate_r1}"
            RAW_R2="${candidate_r2}"
            return 0
        fi
    done
    fail_step "input" 2 "在 ${DATA_DIR} 中未找到样本 ${SAMPLE} 的双端 FASTQ 文件"
}

# 解析工具输出的双端文件（按常见命名模式匹配）
resolve_pair() {
    local prefix="$1" label="$2"
    local pairs=(
        "${prefix}_fwd.fq.gz:${prefix}_rev.fq.gz"
        "${prefix}_fwd.fastq.gz:${prefix}_rev.fastq.gz"
        "${prefix}_fwd.fq:${prefix}_rev.fq"
        "${prefix}_fwd.fastq:${prefix}_rev.fastq"
        "${prefix}_1.fq.gz:${prefix}_2.fq.gz"
        "${prefix}_1.fastq.gz:${prefix}_2.fastq.gz"
        "${prefix}_1.fq:${prefix}_2.fq"
        "${prefix}_1.fastq:${prefix}_2.fastq"
        "${prefix}.1.fq.gz:${prefix}.2.fq.gz"
        "${prefix}.1.fastq.gz:${prefix}.2.fastq.gz"
        "${prefix}.1.fq:${prefix}.2.fq"
        "${prefix}.1.fastq:${prefix}.2.fastq"
        "${prefix}_R1.fq.gz:${prefix}_R2.fq.gz"
        "${prefix}_R1.fastq.gz:${prefix}_R2.fastq.gz"
        "${prefix}_R1.fq:${prefix}_R2.fq"
        "${prefix}_R1.fastq:${prefix}_R2.fastq"
    )

    if [ "${DRY_RUN}" = "true" ]; then
        RESOLVED_R1="${pairs[0]%%:*}"
        RESOLVED_R2="${pairs[0]##*:}"
        return 0
    fi

    for item in "${pairs[@]}"; do
        left="${item%%:*}"
        right="${item##*:}"
        if [ -s "${left}" ] && [ -s "${right}" ]; then
            RESOLVED_R1="${left}"
            RESOLVED_R2="${right}"
            return 0
        fi
    done
    fail_step "${label}" 2 "无法解析双端输出文件: ${prefix}"
}

# 按最小氨基酸长度过滤 Prodigal 输出的蛋白 FASTA
filter_faa_by_length() {
    local input="$1" output="$2" min_aa="$3"
    "${PYTHON_BIN}" - "$input" "$output" "$min_aa" <<'PY'
import sys

src, dst, min_aa = sys.argv[1], sys.argv[2], int(sys.argv[3])
kept = total = 0
name, seq = None, []
with open(src, encoding="utf-8") as handle, open(dst, "w", encoding="utf-8") as out:
    def flush():
        global kept, total
        if name is None:
            return
        total += 1
        joined = "".join(seq)
        if len(joined) >= min_aa:
            kept += 1
            out.write(name + "\n")
            for i in range(0, len(joined), 60):
                out.write(joined[i:i+60] + "\n")
    for raw in handle:
        line = raw.rstrip()
        if line.startswith(">"):
            flush()
            name, seq = line, []
        else:
            seq.append(line)
    flush()
print(f"filtered_proteins kept={kept} total={total} min_aa={min_aa}")
PY
}

# ============================================================
# 预检查：验证必需工具
# ============================================================
if [ "${DRY_RUN}" != "true" ]; then
    require_tool fastp
    require_tool "prinseq++"
    require_tool bowtie2
    require_tool megahit
    require_tool prodigal
fi

# ============================================================
# 查找输入文件
# ============================================================
find_reads
log_msg "样本=${SAMPLE} 配置=${CONFIG_YAML} 每进程线程=${PER_PROC_THREADS} 总线程=${THREADS}"
log_msg "输入 R1=${RAW_R1} R2=${RAW_R2}"

# ============================================================
# 步骤 1: fastp — 质量修剪 + 接头处理
# ============================================================
FASTP_R1="${FASTP_DIR}/${SAMPLE}_fastp_R1.fq.gz"
FASTP_R2="${FASTP_DIR}/${SAMPLE}_fastp_R2.fq.gz"
CMD=(fastp -i "${RAW_R1}" -I "${RAW_R2}" -o "${FASTP_R1}" -O "${FASTP_R2}")

# 接头/修剪参数
[ "${FASTP_DISABLE_ADAPTER}" = "true" ] && CMD+=(--disable_adapter_trimming)
[ "${FASTP_DISABLE_ADAPTER}" != "true" ] && [ "${FASTP_DETECT_ADAPTER}" = "true" ] && CMD+=(--detect_adapter_for_pe)
[ "${FASTP_CUT_FRONT}" = "true" ] && CMD+=(--cut_front --cut_front_window_size "${FASTP_FRONT_WINDOW}" --cut_front_mean_quality "${FASTP_FRONT_QUAL}")
[ "${FASTP_CUT_TAIL}" = "true" ]  && CMD+=(--cut_tail  --cut_tail_window_size  "${FASTP_TAIL_WINDOW}"  --cut_tail_mean_quality  "${FASTP_TAIL_QUAL}")
[ "${FASTP_CUT_RIGHT}" = "true" ] && CMD+=(--cut_right --cut_right_window_size "${FASTP_RIGHT_WINDOW}" --cut_right_mean_quality "${FASTP_RIGHT_QUAL}")

CMD+=(--qualified_quality_phred "${FASTP_QUALIFIED_PHRED}" --unqualified_percent_limit "${FASTP_UNQUALIFIED_PERCENT}")
CMD+=(--length_required "${FASTP_LENGTH}" --n_base_limit "${FASTP_N_LIMIT}")
CMD+=(-h "${FASTP_DIR}/${SAMPLE}_fastp.html" -j "${FASTP_DIR}/${SAMPLE}_fastp.json" -w "${PER_PROC_THREADS}")
append_extra_args "${FASTP_EXTRA_ARGS}"
run_step fastp "${CMD[@]}"

# ============================================================
# 步骤 2: PRINSEQ++ — 低复杂度序列过滤（PDF p23）
# ============================================================
PRINSEQ_PREFIX="${PRINSEQ_DIR}/${SAMPLE}_prinseq"
CMD=(prinseq++
    -fastq "${FASTP_R1}" -fastq2 "${FASTP_R2}"
    -out_name "${PRINSEQ_PREFIX}"
    -lc_entropy "${PRINSEQ_ENTROPY}" -lc_dust "${PRINSEQ_DUST}"
    -threads "${PER_PROC_THREADS}")
append_extra_args "${PRINSEQ_EXTRA_ARGS}"
run_step prinseq "${CMD[@]}"
resolve_pair "${PRINSEQ_PREFIX}_good_out" "prinseq"
PRINSEQ_R1="${RESOLVED_R1}"
PRINSEQ_R2="${RESOLVED_R2}"

# ============================================================
# 步骤 3: bowtie2 — 去除 rRNA（替换 SortMeRNA；PDF p25）
#   将 reads 比对到 rRNA 参考索引，提取未比对（非 rRNA）reads
# ============================================================
CLEAN_PREFIX="${CLEAN_DIR}/${SAMPLE}_clean_non_rRNA"
CLEAN_R1="${CLEAN_PREFIX}_R1.fq.gz"
CLEAN_R2="${CLEAN_PREFIX}_R2.fq.gz"

CMD=(bowtie2
    -x "${BOWTIE2_INDEX}"
    -1 "${PRINSEQ_R1}" -2 "${PRINSEQ_R2}"
    --un-conc-gz "${CLEAN_PREFIX}"
    "--${BOWTIE2_PRESET}"
    -p "${PER_PROC_THREADS}"
    --no-unal
    -S /dev/null)
append_extra_args "${BOWTIE2_EXTRA_ARGS}"
run_step bowtie2 "${CMD[@]}"

# bowtie2 --un-conc-gz 输出为 <prefix>.1 和 <prefix>.2（gzipped），重命名为标准双端命名
if [ "${DRY_RUN}" != "true" ]; then
    if [ -f "${CLEAN_PREFIX}.1" ] && [ -f "${CLEAN_PREFIX}.2" ]; then
        mv "${CLEAN_PREFIX}.1" "${CLEAN_R1}"
        mv "${CLEAN_PREFIX}.2" "${CLEAN_R2}"
    else
        fail_step bowtie2 1 "bowtie2 未产生非 rRNA reads 输出"
    fi
    [ -s "${CLEAN_R1}" ] && [ -s "${CLEAN_R2}" ] || \
        fail_step bowtie2 1 "去 rRNA 后无有效 reads 剩余"
fi

log_msg "bowtie2 去 rRNA 完成: ${CLEAN_R1} / ${CLEAN_R2}"

# ============================================================
# 步骤 4: MEGAHIT — 从头组装（PDF p23-24: 默认参数）
# ============================================================
MEGAHIT_OUT="${MEGAHIT_DIR}/${SAMPLE}"

# 处理已存在的输出目录
if [ -d "${MEGAHIT_OUT}" ]; then
    if [ "${MEGAHIT_OVERWRITE}" = "true" ]; then
        case "${MEGAHIT_OUT}" in
            "${MEGAHIT_DIR}/${SAMPLE}") rm -rf -- "${MEGAHIT_OUT}" ;;
            *) fail_step megahit 2 "拒绝删除非预期的 MEGAHIT 目录: ${MEGAHIT_OUT}" ;;
        esac
    else
        fail_step megahit 2 "MEGAHIT 输出已存在; 设置 megahit.overwrite: true 或手动删除 ${MEGAHIT_OUT}"
    fi
fi

CMD=(megahit
    -1 "${CLEAN_R1}" -2 "${CLEAN_R2}"
    -o "${MEGAHIT_OUT}"
    --min-contig-len "${MEGAHIT_MIN_CONTIG_LEN}"
    --num-cpu-threads "${PER_PROC_THREADS}")
append_extra_args "${MEGAHIT_EXTRA_ARGS}"
run_step megahit "${CMD[@]}"

# 验证组装结果
CONTIGS_FASTA="${MEGAHIT_OUT}/final.contigs.fa"
if [ "${DRY_RUN}" != "true" ]; then
    [ -s "${CONTIGS_FASTA}" ] || fail_step megahit 1 "MEGAHIT 未产生 final.contigs.fa"
fi

# ============================================================
# 可选清理：MEGAHIT 完成后删除前三步中间 fastq
#   开关: config.yaml 设置 cleanup.intermediates: true（默认 false）
#   说明: Prodigal 仅依赖 MEGAHIT 的 final.contigs.fa，删除安全；
#         fastp 的 html/json 质控报告保留
# ============================================================
if [ "${CLEANUP_INTERMEDIATES}" = "true" ] && [ "${DRY_RUN}" != "true" ]; then
    rm -f -- "${FASTP_R1}" "${FASTP_R2}"
    rm -f -- "${PRINSEQ_PREFIX}"*
    rm -f -- "${CLEAN_R1}" "${CLEAN_R2}"
    # 清空路径，避免 JSON outputs 中记录已删除文件
    FASTP_R1=""; FASTP_R2=""
    PRINSEQ_R1=""; PRINSEQ_R2=""
    CLEAN_R1=""; CLEAN_R2=""
    log_msg "中间文件清理完成（cleanup.intermediates=true）"
fi

# ============================================================
# 步骤 5: Prodigal — ORF 预测（PDF p24: ≥200 aa 蛋白保留）
# ============================================================
PROTEIN_FAA="${PRODIGAL_DIR}/${SAMPLE}_final.contigs_proteins.faa"
PROTEIN_FFN="${PRODIGAL_DIR}/${SAMPLE}_final.contigs_genes.ffn"
PRODIGAL_GFF="${PRODIGAL_DIR}/${SAMPLE}_final.contigs.gff"

CMD=(prodigal
    -i "${CONTIGS_FASTA}"
    -a "${PROTEIN_FAA}" -d "${PROTEIN_FFN}" -o "${PRODIGAL_GFF}"
    -f gff -p "${PRODIGAL_MODE}" -q)
append_extra_args "${PRODIGAL_EXTRA_ARGS}"
run_step prodigal "${CMD[@]}"

[ "${DRY_RUN}" = "true" ] || [ -s "${PROTEIN_FAA}" ] || \
    fail_step prodigal 1 "Prodigal 未产生蛋白序列"

# 过滤 ≥200 aa 的蛋白（PDF: 保留 ≥200 氨基酸残基的 ORF）
FILTERED_FAA="${PRODIGAL_DIR}/${SAMPLE}_final.contigs_proteins.min${PRODIGAL_MIN_AA}aa.faa"
if [ "${DRY_RUN}" != "true" ]; then
    filter_faa_by_length "${PROTEIN_FAA}" "${FILTERED_FAA}" "${PRODIGAL_MIN_AA}" \
        >> "${LOG_DIR}/prodigal.log" 2>&1
fi

# ============================================================
# 完成：写入成功 JSON，管道结束
# ============================================================
write_json "success" "" "流程完成" 0
