#!/usr/bin/env bash
# ============================================================
# 4.2 — LucaProt RdRp probability prediction
#
# A batch step: it scores every L-segment protein collected in
# <results_dir>/00_results/l.fasta, which the tree step (4-vmr-tree.py)
# writes. The probability of each protein being an RdRp becomes the third
# line of evidence used by the risk-scoring agent.
#
# Input
#   <results_dir>/00_results/l.fasta
#
# Output
#   <results_dir>/00_results/l.fasta.csv   id + RdRp probability
#   <results_dir>/json/lucaprot.json       structured result
#   <results_dir>/log/lucaprot.log         full command and tool output
#
# Usage
#   bash 4.2-lucaprot.sh -c ../../config.yaml
#   bash 4.2-lucaprot.sh -c ../../config.yaml --dry-run
#   bash 4.2-lucaprot.sh -c ../../config.yaml --overwrite
# ============================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON:-python3}"

usage() {
    cat <<'EOF'
Usage:
  bash 4.2-lucaprot.sh [-c <config.yaml>] [--dry-run] [--overwrite]

Options:
  -c, --config <path>   Pipeline configuration. When omitted, config.yaml is
                        searched upward from this script (four levels).
      --dry-run         Resolve the command and write the log, then stop.
                        Nothing is deleted and no prediction is run.
      --overwrite       Replace an existing l.fasta.csv.
  -h, --help            Show this message.
EOF
}

# ---------------------------------------------------------------- arguments ---
CONFIG_YAML=""
DRY_RUN="false"
OVERWRITE="false"

while [ "$#" -gt 0 ]; do
    case "$1" in
        -c|--config) CONFIG_YAML="${2:-}"; shift 2 ;;
        --dry-run)   DRY_RUN="true"; shift ;;
        --overwrite) OVERWRITE="true"; shift ;;
        -h|--help)   usage; exit 0 ;;
        -*)          echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
        *)           echo "Unexpected argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# ---------------------------------------------- locate the configuration -----
find_default_config() {
    local d="$SCRIPT_DIR"
    local name
    for _ in 1 2 3 4; do
        for name in config.yaml config-cell-paper.yaml cell_submit-paper.yaml; do
            if [ -f "${d}/${name}" ]; then
                printf '%s\n' "${d}/${name}"
                return 0
            fi
        done
        d="$(dirname "${d}")"
    done
    return 1
}

if [ -z "${CONFIG_YAML}" ]; then
    if ! CONFIG_YAML="$(find_default_config)"; then
        echo "No configuration file found. Pass -c <config.yaml>." >&2
        exit 2
    fi
fi
[ -f "${CONFIG_YAML}" ] || { echo "Configuration not found: ${CONFIG_YAML}" >&2; exit 2; }
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || {
    echo "python3 is required to parse the configuration." >&2
    exit 127
}
command -v conda >/dev/null 2>&1 || { echo "conda is not on PATH." >&2; exit 127; }

# ------------------------------------------------ read the configuration ----
CONFIG_ENV="$(SKILL_SCRIPT_DIR="${SCRIPT_DIR}" "${PYTHON_BIN}" - "${CONFIG_YAML}" <<'PY'
import os
import sys


def resolve_code_dir(configured, script_dir):
    """Return (code_dir, tried_paths).

    An explicit lucaprot.code_dir always wins. When it is empty the script
    looks for the LucaProt checkout next to this skill, so a fresh clone works
    with no configuration at all.
    """
    if configured:
        return configured, []

    tried = []
    if script_dir:
        # <repo>/skills/<agent>/scripts -> <repo>/scripts/LucaProt/src
        tried.append(os.path.join(script_dir, "..", "..", "..", "scripts",
                                  "LucaProt", "src"))
        # <repo>/skills/<agent>/scripts -> <repo>/LucaProt/src
        tried.append(os.path.join(script_dir, "..", "..", "..", "LucaProt", "src"))
        # flat layout: <root>/<script>.py -> <root>/scripts/LucaProt/src
        tried.append(os.path.join(script_dir, "..", "scripts", "LucaProt", "src"))
    for cand in tried:
        cand = os.path.abspath(cand)
        if os.path.isfile(os.path.join(cand, "predict_many_samples.py")):
            return cand, tried
    return "", tried


def read_simple_yaml(path):
    """Minimal indentation-aware YAML reader: scalars and one nesting level."""
    cfg = {}
    section = None
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            indent = len(line) - len(line.lstrip(" "))
            if ":" not in line:
                continue
            key, value = line.strip().split(":", 1)
            key = key.strip()
            value = value.split(" #", 1)[0].strip().strip("\"'")
            if indent == 0:
                if value == "":
                    cfg[key] = {}
                    section = key
                else:
                    cfg[key] = value
                    section = None
            elif section is not None:
                if not isinstance(cfg.get(section), dict):
                    cfg[section] = {}
                cfg[section][key] = value
    return cfg


cfg = read_simple_yaml(sys.argv[1])
project_dir = cfg.get("project_dir", "")
if not project_dir:
    sys.stderr.write("project_dir is missing from the configuration\n")
    sys.exit(2)

paths = cfg.get("paths") or {}
results_dir = paths.get("results_dir") or os.path.join(project_dir, "results")

lp = cfg.get("lucaprot") or {}
code_dir, tried = resolve_code_dir(lp.get("code_dir", ""),
                                   os.environ.get("SKILL_SCRIPT_DIR", ""))
if not code_dir:
    sys.stderr.write("lucaprot.code_dir is empty in the configuration and no "
                     "LucaProt checkout was found automatically.\n")
    sys.stderr.write("Clone it first:\n")
    sys.stderr.write("    git clone https://github.com/alibaba/LucaProt.git "
                     "<repo>/scripts/LucaProt\n")
    if tried:
        sys.stderr.write("Looked for predict_many_samples.py under:\n")
        for cand in tried:
            sys.stderr.write("    %s\n" % os.path.abspath(cand))
    sys.exit(2)

python_bin = lp.get("python", "")
emb_dir = lp.get("emb_dir", "") or os.path.join(code_dir, "..", "emb")

values = {
    "PROJECT_DIR": project_dir,
    "CONDA_ENV": lp.get("conda_env") or cfg.get("conda_env", "lpy-hanta-v2"),
    "RESULTS_DIR": results_dir,
    "LP_CODE_DIR": os.path.abspath(code_dir),
    "LP_PYTHON": python_bin or "python3",
    "LP_EMB_DIR": os.path.abspath(emb_dir),
    "LP_GPU": lp.get("gpu_id", "0"),
    "LP_THRESH": lp.get("threshold", "0.5"),
    "LP_TRUNC": lp.get("truncation_seq_length", "4096"),
    "LP_DATASET": lp.get("dataset_name", "rdrp_40_extend"),
    "LP_DSTYPE": lp.get("dataset_type", "protein"),
    "LP_TASK": lp.get("task_type", "binary_class"),
    "LP_MODEL": lp.get("model_type", "sefn"),
    "LP_TIME": lp.get("time_str", "20230201140320"),
    "LP_STEP": lp.get("step", "100000"),
    "LP_PRINT_N": lp.get("print_per_number", "10"),
    "LP_EXTRA": lp.get("extra_args", ""),
}
for key, value in values.items():
    print("%s=%s" % (key, "'" + str(value).replace("'", "'\\''") + "'"))
PY
)" || { echo "Failed to parse configuration: ${CONFIG_YAML}" >&2; exit 2; }
eval "${CONFIG_ENV}"

[ -d "${LP_CODE_DIR}" ] || { echo "lucaprot.code_dir not found: ${LP_CODE_DIR}" >&2; exit 2; }
[ -f "${LP_CODE_DIR}/predict_many_samples.py" ] || {
    echo "predict_many_samples.py not found under ${LP_CODE_DIR}" >&2
    exit 2
}

# ------------------------------------------------------------- directories ---
OUT_DIR="${RESULTS_DIR}/00_results"
FASTA="${OUT_DIR}/l.fasta"
CSV="${OUT_DIR}/l.fasta.csv"
JSON_DIR="${RESULTS_DIR}/json"
JSON_OUT="${JSON_DIR}/lucaprot.json"
LOG_DIR="${RESULTS_DIR}/log"
LOG_FILE="${LOG_DIR}/lucaprot.log"

mkdir -p "${OUT_DIR}" "${JSON_DIR}" "${LOG_DIR}"

# --------------------------------------------------------------- preflight ---
write_fail_json() {
    local message="$1" code="$2"
    MESSAGE="${message}" CODE="${code}" JSON_OUT="${JSON_OUT}" \
    "${PYTHON_BIN}" - <<'PY'
import json
import os

payload = {
    "lucaprot": {
        "status": "failed",
        "message": os.environ["MESSAGE"],
        "exit_code": int(os.environ["CODE"]),
    }
}
with open(os.environ["JSON_OUT"], "w", encoding="utf-8") as fh:
    fh.write(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n")
PY
}

if [ ! -s "${FASTA}" ]; then
    write_fail_json "input not found or empty: ${FASTA}" 2
    echo "Input not found or empty: ${FASTA}" >&2
    echo "Run the tree step (4-vmr-tree.py) first; it writes l.fasta." >&2
    exit 2
fi

if [ -e "${CSV}" ] && [ "${OVERWRITE}" != "true" ]; then
    write_fail_json "output already exists, pass --overwrite: ${CSV}" 2
    echo "Output already exists, pass --overwrite: ${CSV}" >&2
    exit 2
fi

# ---------------------------------------------------------------- assemble ---
CMD=(conda run -n "${CONDA_ENV}" "${LP_PYTHON}" predict_many_samples.py
     --fasta_file "${FASTA}"
     --save_file "${CSV}"
     --emb_dir "${LP_EMB_DIR}"
     --truncation_seq_length "${LP_TRUNC}"
     --dataset_name "${LP_DATASET}"
     --dataset_type "${LP_DSTYPE}"
     --task_type "${LP_TASK}"
     --model_type "${LP_MODEL}"
     --time_str "${LP_TIME}"
     --step "${LP_STEP}"
     --threshold "${LP_THRESH}"
     --print_per_number "${LP_PRINT_N}"
     --gpu_id "${LP_GPU}")
if [ -n "${LP_EXTRA}" ]; then
    read -r -a _EXTRA <<< "${LP_EXTRA}"
    CMD+=("${_EXTRA[@]}")
fi

{
    printf '[%s] CMD:' "$(date '+%F %T')"
    printf ' %q' "${CMD[@]}"
    printf '\n'
    printf '[%s] WORKDIR: %s\n' "$(date '+%F %T')" "${LP_CODE_DIR}"
} > "${LOG_FILE}"

# -------------------------------------------------------------------- run ----
if [ "${DRY_RUN}" = "true" ]; then
    printf '[%s] DRY_RUN: prediction skipped\n' "$(date '+%F %T')" >> "${LOG_FILE}"
    STATUS="success"
    MESSAGE="dry-run: prediction skipped"
    CODE=0
else
    set +e
    ( cd "${LP_CODE_DIR}" && "${CMD[@]}" ) >> "${LOG_FILE}" 2>&1
    CODE=$?
    set -e
    if [ "${CODE}" -eq 0 ] && [ -s "${CSV}" ]; then
        STATUS="success"
        MESSAGE="success"
    elif [ "${CODE}" -eq 0 ]; then
        STATUS="failed"
        MESSAGE="prediction produced no output file: ${CSV}"
        CODE=1
    else
        STATUS="failed"
        MESSAGE="LucaProt prediction failed, see ${LOG_FILE}"
    fi
fi

# ------------------------------------------------------------ result JSON ----
STATUS="${STATUS}" MESSAGE="${MESSAGE}" CODE="${CODE}" \
FASTA="${FASTA}" CSV="${CSV}" LOG_FILE="${LOG_FILE}" JSON_OUT="${JSON_OUT}" \
LP_THRESH="${LP_THRESH}" LP_GPU="${LP_GPU}" CONDA_ENV="${CONDA_ENV}" \
"${PYTHON_BIN}" - <<'PY'
import json
import os

payload = {
    "lucaprot": {
        "status": os.environ["STATUS"],
        "message": os.environ["MESSAGE"],
        "exit_code": int(os.environ["CODE"]),
        "query": os.environ["FASTA"],
        "outputs": {"predictions": os.environ["CSV"]},
        "logs": {"predict": os.environ["LOG_FILE"]},
        "params": {
            "threshold": float(os.environ["LP_THRESH"]),
            "gpu_id": os.environ["LP_GPU"],
            "conda_env": os.environ["CONDA_ENV"],
        },
    }
}
with open(os.environ["JSON_OUT"], "w", encoding="utf-8") as fh:
    fh.write(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n")
PY

if [ "${STATUS}" != "success" ]; then
    echo "${MESSAGE}" >&2
    exit "${CODE}"
fi

# ------------------------------------------------------------ prediction CB ---
printf '[%s] DONE: %s\n' "$(date '+%F %T')" "${CSV}" >> "${LOG_FILE}"
exit 0
