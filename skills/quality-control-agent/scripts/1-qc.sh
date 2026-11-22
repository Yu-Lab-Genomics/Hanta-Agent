#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# Hantavirus single-sample pre-processing and assembly pipeline
# Steps: fastp -> PRINSEQ++ -> bowtie2 (rRNA removal) -> MEGAHIT -> Prodigal
#
# Usage:
#   bash scripts/qc.sh -s <SAMPLE> [-c config.yaml] [--dry-run]
#   bash scripts/qc.sh -c config.yaml [--dry-run]          (auto-discover every sample in data_dir)
#
# Output:
#   Written to <project_dir>/results by default; override with paths.results_dir in config.yaml
#   <results_dir>/json/<SAMPLE>.qc.json   - structured result, consumed by downstream skills
#   <results_dir>/log/<SAMPLE>/           - one log file per step
#   Silent: nothing is printed to the screen
#
# Optional cleanup:
#   With cleanup.intermediates: true, the intermediate FASTQs of the first
#   three steps are deleted once MEGAHIT finishes (fastp html/json QC reports are kept)
# ============================================================

# ---------- Path initialisation ----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Default configuration: walk up from the script directory looking for config.yaml, so no absolute path is hard-coded across machines
# Layout skills_new/<skill>/scripts/1-qc.sh -> yaml lives at skills_new/config.yaml (two levels up)
# An explicit -c still wins
DEFAULT_CONFIG="${SCRIPT_DIR}/config.yaml"
for _dir in "${SCRIPT_DIR}" "$(dirname "${SCRIPT_DIR}")" "$(dirname "$(dirname "${SCRIPT_DIR}")")"; do
    if [ -f "${_dir}/config.yaml" ]; then
        DEFAULT_CONFIG="${_dir}/config.yaml"
        break
    fi
done
# Backwards compatibility: if config.yaml is absent, fall back to cell_submit-paper.yaml beside the script
if [ ! -f "${DEFAULT_CONFIG}" ] && [ -f "${SCRIPT_DIR}/cell_submit-paper.yaml" ]; then
    DEFAULT_CONFIG="${SCRIPT_DIR}/cell_submit-paper.yaml"
fi
PYTHON_BIN="${PYTHON:-python3}"

SAMPLE=""
CONFIG_YAML="${CONFIG_YAML:-${DEFAULT_CONFIG}}"
DRY_RUN="false"
COMPLETED_STEPS=""

# ---------- Help ----------
usage() {
    cat <<'EOF'
Usage:
    bash scripts/qc.sh -s <sample> [-c <config>] [--dry-run]
    bash scripts/qc.sh <sample>
    bash scripts/qc.sh -c <config> [--dry-run]          (auto-discover every sample in data_dir)
EOF
}

# ---------- JSON error output (written to a file, never printed) ----------
plain_json_error() {
    local sample="${1:-unknown}" message="$2" code="${3:-1}"
    local err_json="${JSON_DIR:-/tmp}/${sample}.qc.json"
    mkdir -p "$(dirname "${err_json}")" 2>/dev/null || true
    printf '{"%s":{"sample":"%s","status":"failed","message":"%s","exit_code":%s}}\n' \
        "${sample}" "${sample}" "${message//\"/\\\"}" "${code}" > "${err_json}"
    exit "${code}"
}

# ---------- Argument parsing ----------
while [ "$#" -gt 0 ]; do
    case "$1" in
        -s|--sample) SAMPLE="${2:-}"; shift 2 ;;
        -c|--config) CONFIG_YAML="${2:-}"; shift 2 ;;
        --dry-run)   DRY_RUN="true"; shift ;;
        -h|--help)   usage; exit 0 ;;
        -*)
            usage >&2
            plain_json_error "${SAMPLE:-unknown}" "unknown option: $1" 2
            ;;
        *)
            if [ -z "${SAMPLE}" ]; then SAMPLE="$1"; shift
            else usage >&2; plain_json_error "${SAMPLE}" "unexpected argument: $1" 2
            fi
            ;;
    esac
done

# ---------- Preflight checks ----------
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || plain_json_error "${SAMPLE:-unknown}" "python3 is required to parse config.yaml" 127
[ -f "${CONFIG_YAML}" ] || plain_json_error "${SAMPLE:-unknown}" "config file not found: ${CONFIG_YAML}" 2

# ============================================================
# Configuration loading: read every parameter from config.yaml and export it as a shell variable
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
    if not isinstance(data, dict):
        return default
    value = data.get(key, "")
    # An explicitly empty value falls back to the default. Without this a key
    # left blank in config.example.yaml, such as ref_dir, would emit an empty
    # string and bowtie2 would silently search the current directory.
    return value if value != "" else default

def emit(name, value):
    print(f"{name}={shlex.quote(str(value))}")

cfg = read_simple_yaml(cfg_path)
project_dir = cfg.get("project_dir")
if not project_dir:
    raise SystemExit("ERROR: project_dir must be set in config.yaml")

paths = cfg.get("paths", {}) if isinstance(cfg.get("paths"), dict) else {}
results_dir = paths.get("results_dir") or os.path.join(project_dir, "results")

# --- Directory paths ---
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

# --- Input discovery ---
emit("R1_SUFFIXES", get("input", "r1_suffixes", "_r1.fq.gz,_R1.fq.gz,_1.fastq.gz"))
emit("R2_SUFFIXES", get("input", "r2_suffixes", "_r2.fq.gz,_R2.fq.gz,_2.fastq.gz"))

# --- Threads ---
emit("THREADS_CFG", cfg.get("threads", "auto"))
# --- Conda environment ---
emit("CONDA_ENV", cfg.get("conda_env", "lpy-hanta-v2"))

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

# --- bowtie2 (replaces SortMeRNA) ---
bt2 = lambda k, d: get("bowtie2", k, d)
emit("BOWTIE2_REF_DIR", bt2("ref_dir", os.path.join(project_dir, "DB", "sortmerna")))
emit("BOWTIE2_INDEX", os.path.join(bt2("ref_dir", os.path.join(project_dir, "DB", "sortmerna")), bt2("index_prefix", "rrna_index")))
emit("BOWTIE2_PRESET", bt2("preset", "very-sensitive"))
emit("BOWTIE2_AUTO_FETCH", bt2("auto_fetch", "true"))
emit("BOWTIE2_FETCH_URL", bt2("fetch_url", "https://github.com/Yu-Lab-Genomics/Hanta-Agent/releases/download/sortmerna-v1.0.0/rrna_index.tar.gz"))
emit("BOWTIE2_FETCH_SHA256", bt2("fetch_sha256", ""))
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

# --- Intermediate cleanup (delete the first three steps' FASTQs after MEGAHIT; off by default) ---
cl = lambda k, d: get("cleanup", k, d)
emit("CLEANUP_INTERMEDIATES", cl("intermediates", "false"))
PY
}

CONFIG_ENV="$(load_config)" || plain_json_error "${SAMPLE:-unknown}" "failed to parse config.yaml: ${CONFIG_YAML}" 2
eval "${CONFIG_ENV}"

# ============================================================
# Thread and parallelism calculation (moved up here so the multi-sample scheduler can use it)
# ============================================================

# Auto-detect available threads (read threads from config.yaml; 'auto' probes the machine)
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

# Threads per process: fastp caps at 16 per process, the other tools are aligned with it
# so the total core count never exceeds the configured value when samples run in parallel
PER_PROC_THREADS=16
# Parallel samples = total threads / threads per process (floor, at least 1)
# Guards against a total that is not a multiple of 16: the remainder cores stay idle, never oversubscribed
if [ "${THREADS}" -gt "${PER_PROC_THREADS}" ]; then
    PARALLEL=$(( THREADS / PER_PROC_THREADS ))
else
    PARALLEL=1
    PER_PROC_THREADS="${THREADS}"
fi
[ "${PARALLEL}" -ge 1 ] || PARALLEL=1

# ============================================================
# Auto-discover samples (scan data_dir when -s is not given)
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
    # Deduplicate
    SAMPLES_FOUND="$(printf '%s\n' ${SAMPLES_FOUND} | sort -u | tr '\n' ' ')"
    SAMPLES_FOUND="${SAMPLES_FOUND% }"

    if [ -z "${SAMPLES_FOUND}" ]; then
        usage >&2
        plain_json_error "unknown" "no sample could be auto-discovered in ${DATA_DIR} and -s was not given (suffixes: ${R1_SUFFIXES})" 2
    fi

    printf '[%s] auto-discovered %s sample(s), %s in parallel (%s threads per sample per tool)\n' \
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
# Utility functions
# ============================================================

RESULT_JSON="${JSON_DIR}/${SAMPLE}.qc.json"

# Create every output directory in one pass
mkdir -p "${FASTP_DIR}" "${PRINSEQ_DIR}" "${CLEAN_DIR}" \
         "${MEGAHIT_DIR}" "${PRODIGAL_DIR}" "${LOG_DIR}" "${JSON_DIR}"

# ---------- JSON result writing (file only, never printed) ----------
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
# Never printed to stdout
PY
}

# Called on failure: write an error JSON and exit
fail_step() {
    local step="$1" code="$2" message="$3"
    # Also report to stderr. The JSON below records the reason for tooling, but
    # a person watching the terminal only sees this line.
    printf '[%s] FAILED %s: %s\n' "$(date '+%F %T')" "${step}" "${message}" >&2
    write_json "failed" "${step}" "${message}" "${code}"
    exit "${code}"
}

# Append to the pipeline log
log_msg() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" >> "${LOG_DIR}/pipeline.log"
}

# Check that a required tool exists
require_tool() {
    local tool="$1"
    command -v "${tool}" >/dev/null 2>&1 || fail_step "preflight" 127 "required tool not found: ${tool}"
}

# ---------------------------------------------------------------
# rRNA index: fetch it on demand
#   The six bowtie2 files total about 1.2 GB and are published as a release
#   asset rather than tracked in git, so a fresh clone has none of them.
#   Downloading here means the pipeline can start straight after cloning, with
#   no manual preparation step.
# ---------------------------------------------------------------
RRNA_INDEX_FILES=(
    rrna_index.1.bt2
    rrna_index.2.bt2
    rrna_index.3.bt2
    rrna_index.4.bt2
    rrna_index.rev.1.bt2
    rrna_index.rev.2.bt2
)

# How many of the six files are absent under the configured directory
rrna_index_missing_count() {
    local f n=0
    for f in "${RRNA_INDEX_FILES[@]}"; do
        [ -s "${BOWTIE2_REF_DIR}/${f}" ] || n=$((n + 1))
    done
    printf '%s\n' "${n}"
}

# Printed by every failure path, so the wording cannot drift between them
rrna_index_hint() {
    printf '  Prepare the index by hand with:\n' >&2
    printf '    bash skills/database-update-agent/scripts/fetch_rrna_index.sh --db-dir %s/DB\n' "${PROJECT_DIR}" >&2
    printf '  Or set bowtie2.auto_fetch: true in the configuration.\n' >&2
}

fetch_rrna_index() {
    local url="${BOWTIE2_FETCH_URL}"
    local want_sha="${BOWTIE2_FETCH_SHA256:-}"
    local dir="${BOWTIE2_REF_DIR}"
    local name archive tmp f src got

    name="$(basename "${url}")"
    archive="${dir}/${name}"

    mkdir -p "${dir}" || fail_step "preflight" 1 "cannot create the rRNA index directory: ${dir}"

    printf '[%s] rRNA index incomplete under %s\n' "$(date '+%F %T')" "${dir}"
    printf '[%s] downloading %s (about 1.2 GB, once only)\n' "$(date '+%F %T')" "${url}"
    log_msg "rRNA index auto-fetch: ${url}"

    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 --retry-delay 5 -o "${archive}.partial" "${url}" || {
            rm -f "${archive}.partial"
            rrna_index_hint
            fail_step "preflight" 1 "rRNA index download failed: ${url}"
        }
    elif command -v wget >/dev/null 2>&1; then
        wget --tries=3 --waitretry=5 -O "${archive}.partial" "${url}" || {
            rm -f "${archive}.partial"
            rrna_index_hint
            fail_step "preflight" 1 "rRNA index download failed: ${url}"
        }
    else
        rrna_index_hint
        fail_step "preflight" 1 "neither curl nor wget is available; cannot fetch the rRNA index"
    fi
    mv -f "${archive}.partial" "${archive}"

    if [ -n "${want_sha}" ]; then
        got="$(sha256sum "${archive}" | awk '{print $1}')"
        if [ "${got}" != "${want_sha}" ]; then
            printf '  archive kept at %s for inspection\n' "${archive}" >&2
            fail_step "preflight" 1 "rRNA index checksum mismatch: expected ${want_sha}, got ${got}"
        fi
        printf '[%s] rRNA index checksum OK\n' "$(date '+%F %T')"
    fi

    tmp="$(mktemp -d "${dir}/.extract-XXXXXX")" \
        || fail_step "preflight" 1 "cannot create a temporary directory under ${dir}"
    if ! tar -xzf "${archive}" -C "${tmp}"; then
        fail_step "preflight" 1 "rRNA index extraction failed; leftovers kept in ${tmp}"
    fi

    # Locate each file by name, so an archive that wraps everything in a
    # directory behaves the same as one that does not.
    for f in "${RRNA_INDEX_FILES[@]}"; do
        src="$(find "${tmp}" -type f -name "${f}" -print -quit)"
        if [ -z "${src}" ]; then
            fail_step "preflight" 1 "rRNA index archive does not contain ${f}; extracted tree kept in ${tmp}"
        fi
        mv -f "${src}" "${dir}/${f}"
    done
    rm -rf "${tmp}"
    rm -f "${archive}"

    printf '[%s] rRNA index installed into %s\n' "$(date '+%F %T')" "${dir}"
    log_msg "rRNA index installed into ${dir}"
}

ensure_rrna_index() {
    local missing
    missing="$(rrna_index_missing_count)"
    [ "${missing}" -eq 0 ] && return 0

    if [ "${DRY_RUN}" = "true" ]; then
        printf '[dry-run] rRNA index incomplete: %s of %s files missing from %s\n' \
            "${missing}" "${#RRNA_INDEX_FILES[@]}" "${BOWTIE2_REF_DIR}"
        printf '[dry-run] would download %s and extract it there\n' "${BOWTIE2_FETCH_URL}"
        return 0
    fi

    if [ "${BOWTIE2_AUTO_FETCH}" != "true" ]; then
        printf '[%s] rRNA index incomplete: %s of %s files missing from %s\n' \
            "$(date '+%F %T')" "${missing}" "${#RRNA_INDEX_FILES[@]}" "${BOWTIE2_REF_DIR}" >&2
        rrna_index_hint
        fail_step "preflight" 1 "rRNA index incomplete and bowtie2.auto_fetch is not true"
    fi

    fetch_rrna_index

    [ "$(rrna_index_missing_count)" -eq 0 ] \
        || fail_step "preflight" 1 "rRNA index still incomplete after fetching"
}

# Append extra arguments (split on whitespace)
append_extra_args() {
    local extra="$1"
    [ -n "${extra}" ] || return 0
    read -r -a EXTRA_ARGV <<< "${extra}"
    CMD+=("${EXTRA_ARGV[@]}")
}

# Run one step: log the command, run it, then check the exit code
run_step() {
    local step="$1"; shift
    local log_file="${LOG_DIR}/${step}.log"
    log_msg "START ${step}"

    # Record the command in the log
    {
        printf '[%s] CMD:' "$(date '+%F %T')"
        printf ' %q' "$@"
        printf '\n'
    } > "${log_file}"

    if [ "${DRY_RUN}" = "true" ]; then
        printf '[%s] DRY_RUN: skipping execution\n' "$(date '+%F %T')" >> "${log_file}"
    else
        set +e
        "$@" >> "${log_file}" 2>&1
        local code=$?
        set -e
        if [ "${code}" -ne 0 ]; then
            log_msg "FAIL ${step} exit=${code}; log: ${log_file}"
            fail_step "${step}" "${code}" "${step} failed; see ${log_file}"
        fi
    fi

    COMPLETED_STEPS="${COMPLETED_STEPS:+${COMPLETED_STEPS},}${step}"
    log_msg "DONE ${step}"
}

# Locate the paired-end FASTQ files of a sample inside the data directory
find_reads() {
    local r1_suffix r2_suffix candidate_r1 candidate_r2
    IFS=',' read -r -a R1_ARR <<< "${R1_SUFFIXES}"
    IFS=',' read -r -a R2_ARR <<< "${R2_SUFFIXES}"
    [ "${#R1_ARR[@]}" -eq "${#R2_ARR[@]}" ] || \
        fail_step "input" 2 "r1_suffixes and r2_suffixes must have the same number of entries"

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
    fail_step "input" 2 "no paired-end FASTQ files found for sample ${SAMPLE} in ${DATA_DIR}"
}

# Parse the paired-end outputs of a tool (matched against the common naming patterns)
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
    fail_step "${label}" 2 "cannot parse paired-end output files: ${prefix}"
}

# Filter the Prodigal protein FASTA by minimum amino-acid length
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
# Preflight: verify the required tools
# ============================================================
if [ "${DRY_RUN}" != "true" ]; then
    command -v conda >/dev/null 2>&1 || fail_step "preflight" 127 "conda command not found"
    conda run -n "${CONDA_ENV}" true >/dev/null 2>&1 || fail_step "preflight" 127 "conda environment ${CONDA_ENV} is not usable"
fi

# The rRNA index is the only database component not tracked in git, so it is
# fetched here on demand. Runs in dry-run mode too, where it only reports.
ensure_rrna_index

# ============================================================
# Locate the input files
# ============================================================
find_reads
log_msg "sample=${SAMPLE} config=${CONFIG_YAML} threads_per_process=${PER_PROC_THREADS} total_threads=${THREADS}"
log_msg "input R1=${RAW_R1} R2=${RAW_R2}"

# ============================================================
# Step 1: fastp - quality trimming and adapter removal
# ============================================================
FASTP_R1="${FASTP_DIR}/${SAMPLE}_fastp_R1.fq.gz"
FASTP_R2="${FASTP_DIR}/${SAMPLE}_fastp_R2.fq.gz"
CMD=(conda run -n "${CONDA_ENV}" fastp -i "${RAW_R1}" -I "${RAW_R2}" -o "${FASTP_R1}" -O "${FASTP_R2}")

# Adapter / trimming parameters
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
# Step 2: PRINSEQ++ - low-complexity filtering (PDF p23)
# ============================================================
PRINSEQ_PREFIX="${PRINSEQ_DIR}/${SAMPLE}_prinseq"
CMD=(conda run -n "${CONDA_ENV}" prinseq++
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
# Step 3: bowtie2 - rRNA removal (replaces SortMeRNA; PDF p25)
#   Align reads to the rRNA index and keep the reads that did not align (non-rRNA)
# ============================================================
CLEAN_PREFIX="${CLEAN_DIR}/${SAMPLE}_clean_non_rRNA"
CLEAN_R1="${CLEAN_PREFIX}_R1.fq.gz"
CLEAN_R2="${CLEAN_PREFIX}_R2.fq.gz"

CMD=(conda run -n "${CONDA_ENV}" bowtie2
    -x "${BOWTIE2_INDEX}"
    -1 "${PRINSEQ_R1}" -2 "${PRINSEQ_R2}"
    --un-conc-gz "${CLEAN_PREFIX}"
    "--${BOWTIE2_PRESET}"
    -p "${PER_PROC_THREADS}"
    --no-unal
    -S /dev/null)
append_extra_args "${BOWTIE2_EXTRA_ARGS}"
run_step bowtie2 "${CMD[@]}"

# bowtie2 --un-conc-gz writes <prefix>.1 and <prefix>.2 (gzipped); rename them to the standard paired-end names
if [ "${DRY_RUN}" != "true" ]; then
    if [ -f "${CLEAN_PREFIX}.1" ] && [ -f "${CLEAN_PREFIX}.2" ]; then
        mv "${CLEAN_PREFIX}.1" "${CLEAN_R1}"
        mv "${CLEAN_PREFIX}.2" "${CLEAN_R2}"
    else
        fail_step bowtie2 1 "bowtie2 produced no non-rRNA read output"
    fi
    [ -s "${CLEAN_R1}" ] && [ -s "${CLEAN_R2}" ] || \
        fail_step bowtie2 1 "no usable reads remain after rRNA removal"
fi

log_msg "bowtie2 rRNA removal done: ${CLEAN_R1} / ${CLEAN_R2}"

# ============================================================
# Step 4: MEGAHIT - de novo assembly (PDF p23-24: default parameters)
# ============================================================
MEGAHIT_OUT="${MEGAHIT_DIR}/${SAMPLE}"

# Handle an existing output directory
# Note: deletion happens only on a real run. --dry-run must keep the disk
#       read-only, otherwise a single rehearsal permanently destroys existing assemblies.
if [ -d "${MEGAHIT_OUT}" ]; then
    if [ "${MEGAHIT_OVERWRITE}" = "true" ]; then
        if [ "${DRY_RUN}" = "true" ]; then
            log_msg "DRY-RUN: keeping the existing MEGAHIT directory (a real run would overwrite it): ${MEGAHIT_OUT}"
        else
            case "${MEGAHIT_OUT}" in
                "${MEGAHIT_DIR}/${SAMPLE}") rm -rf -- "${MEGAHIT_OUT}" ;;
                *) fail_step megahit 2 "refusing to delete an unexpected MEGAHIT directory: ${MEGAHIT_OUT}" ;;
            esac
        fi
    else
        fail_step megahit 2 "MEGAHIT output already exists; set megahit.overwrite: true or delete ${MEGAHIT_OUT} manually"
    fi
fi

CMD=(conda run -n "${CONDA_ENV}" megahit
    -1 "${CLEAN_R1}" -2 "${CLEAN_R2}"
    -o "${MEGAHIT_OUT}"
    --min-contig-len "${MEGAHIT_MIN_CONTIG_LEN}"
    --num-cpu-threads "${PER_PROC_THREADS}")
append_extra_args "${MEGAHIT_EXTRA_ARGS}"
run_step megahit "${CMD[@]}"

# Verify the assembly
CONTIGS_FASTA="${MEGAHIT_OUT}/final.contigs.fa"
if [ "${DRY_RUN}" != "true" ]; then
    [ -s "${CONTIGS_FASTA}" ] || fail_step megahit 1 "MEGAHIT produced no final.contigs.fa"
fi

# ============================================================
# Optional cleanup: delete the first three steps' intermediate FASTQs after MEGAHIT
#   Switch: cleanup.intermediates: true in config.yaml (default false)
#   Rationale: Prodigal only needs MEGAHIT's final.contigs.fa, so deletion is safe;
#              the fastp html/json QC reports are kept
# ============================================================
if [ "${CLEANUP_INTERMEDIATES}" = "true" ] && [ "${DRY_RUN}" != "true" ]; then
    rm -f -- "${FASTP_R1}" "${FASTP_R2}"
    rm -f -- "${PRINSEQ_PREFIX}"*
    rm -f -- "${CLEAN_R1}" "${CLEAN_R2}"
    # Clear the paths so the JSON outputs do not point at deleted files
    FASTP_R1=""; FASTP_R2=""
    PRINSEQ_R1=""; PRINSEQ_R2=""
    CLEAN_R1=""; CLEAN_R2=""
    log_msg "intermediate cleanup done (cleanup.intermediates=true)"
fi

# ============================================================
# Step 5: Prodigal - ORF prediction (PDF p24: keep proteins of 200 aa or more)
# ============================================================
PROTEIN_FAA="${PRODIGAL_DIR}/${SAMPLE}_final.contigs_proteins.faa"
PROTEIN_FFN="${PRODIGAL_DIR}/${SAMPLE}_final.contigs_genes.ffn"
PRODIGAL_GFF="${PRODIGAL_DIR}/${SAMPLE}_final.contigs.gff"

CMD=(conda run -n "${CONDA_ENV}" prodigal
    -i "${CONTIGS_FASTA}"
    -a "${PROTEIN_FAA}" -d "${PROTEIN_FFN}" -o "${PRODIGAL_GFF}"
    -f gff -p "${PRODIGAL_MODE}" -q)
append_extra_args "${PRODIGAL_EXTRA_ARGS}"
run_step prodigal "${CMD[@]}"

[ "${DRY_RUN}" = "true" ] || [ -s "${PROTEIN_FAA}" ] || \
    fail_step prodigal 1 "Prodigal produced no protein sequences"

# Keep proteins of 200 aa or more (PDF: retain ORFs of at least 200 residues)
FILTERED_FAA="${PRODIGAL_DIR}/${SAMPLE}_final.contigs_proteins.min${PRODIGAL_MIN_AA}aa.faa"
if [ "${DRY_RUN}" != "true" ]; then
    filter_faa_by_length "${PROTEIN_FAA}" "${FILTERED_FAA}" "${PRODIGAL_MIN_AA}" \
        >> "${LOG_DIR}/prodigal.log" 2>&1
fi

# ============================================================
# Done: write the success JSON and end the pipeline
# ============================================================
write_json "success" "" "pipeline complete" 0
