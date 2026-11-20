#!/usr/bin/env bash
# Fetch the bowtie2 rRNA index used by 1-qc.sh for rRNA depletion.
#
# The index is ~1.2 GB and cannot live in the git repository, so it is published
# as release assets. This script downloads the six .bt2 files, verifies each one
# against its SHA256 checksum, and installs them into <db>/sortmerna/.
#
# Only the .bt2 files are needed. bowtie2 is invoked with -x <prefix>, which
# reads the index and never touches the source FASTA. rrna_merged.fasta is
# therefore published optionally; pass --with-fasta if you also want it for
# rebuilding the index later.
#
# Usage
#   bash fetch_rrna_index.sh --url <release-base-url>
#   bash fetch_rrna_index.sh --url <release-base-url> --db-dir /data/DB
#   bash fetch_rrna_index.sh --db-dir /data/DB --check        # verify only
#   bash fetch_rrna_index.sh --db-dir /data/DB --make-checksums
#
# Exit codes
#   0  success, or the index is already installed and intact
#   1  a download, checksum, or installation step failed
#   2  bad arguments or an unresolved database directory
#
# Every failure is reported on stderr with the path involved. Nothing is deleted:
# a partially downloaded file is left in place under a .partial name so a rerun
# can resume or be inspected.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_NAMES="config.yaml config-cell-paper.yaml cell_submit-paper.yaml"

INDEX_FILES=(
    rrna_index.1.bt2
    rrna_index.2.bt2
    rrna_index.3.bt2
    rrna_index.4.bt2
    rrna_index.rev.1.bt2
    rrna_index.rev.2.bt2
)

DRY_RUN="false"
CHECK_ONLY="false"
MAKE_CHECKSUMS="false"
WITH_FASTA="false"
FORCE="false"
BASE_URL=""
DB_DIR=""
CONFIG_YAML=""

usage() {
    cat <<'EOF'
Usage:
  bash fetch_rrna_index.sh --url <release-base-url> [options]

Options:
  --url <base-url>     Base URL holding the release assets, e.g.
                       https://github.com/OWNER/REPO/releases/download/db-v1.0.0
  --db-dir <path>      Database directory. Defaults to <project_dir>/DB, read
                       from a configuration file found upward from this script.
  -c, --config <path>  Explicit configuration file.
  --check              Verify an already installed index and exit. Writes nothing.
  --make-checksums     Write SHA256SUMS for the files present in <db>/sortmerna.
                       Run this before uploading the assets.
  --with-fasta         Also download rrna_merged.fasta (only needed to rebuild
                       the index; the pipeline never reads it).
  --force              Re-download even when the index looks intact.
  --dry-run            Print what would happen and exit.
  -h, --help           Show this message.
EOF
}

log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
fail() { printf '[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2; exit 1; }

# ---------------------------------------------------------------- arguments ---
while [ "$#" -gt 0 ]; do
    case "$1" in
        --url)            BASE_URL="${2:-}"; shift 2 ;;
        --db-dir)         DB_DIR="${2:-}"; shift 2 ;;
        -c|--config)      CONFIG_YAML="${2:-}"; shift 2 ;;
        --check)          CHECK_ONLY="true"; shift ;;
        --make-checksums) MAKE_CHECKSUMS="true"; shift ;;
        --with-fasta)     WITH_FASTA="true"; shift ;;
        --force)          FORCE="true"; shift ;;
        --dry-run)        DRY_RUN="true"; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)                printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

# ---------------------------------------------- locate the configuration -----
find_default_config() {
    local d="$SCRIPT_DIR" name
    for _ in 1 2 3 4; do
        for name in $CONFIG_NAMES; do
            if [ -f "${d}/${name}" ]; then
                printf '%s\n' "${d}/${name}"
                return 0
            fi
        done
        d="$(dirname "$d")"
    done
    return 1
}

if [ -z "${DB_DIR}" ]; then
    if [ -z "${CONFIG_YAML}" ]; then
        if ! CONFIG_YAML="$(find_default_config)"; then
            fail "no configuration file found. Pass --db-dir or -c <config.yaml>."
        fi
    fi
    [ -f "${CONFIG_YAML}" ] || fail "config file not found: ${CONFIG_YAML}"
    PROJECT_DIR="$(awk '
        /^[[:space:]]*#/ { next }
        /^project_dir[[:space:]]*:/ {
            sub(/^project_dir[[:space:]]*:[[:space:]]*/, "")
            sub(/[[:space:]]+#.*$/, "")
            gsub(/^["\x27]|["\x27]$/, "")
            print
            exit
        }' "$CONFIG_YAML")"
    [ -n "${PROJECT_DIR}" ] || fail "project_dir is missing from ${CONFIG_YAML}"
    DB_DIR="${PROJECT_DIR}/DB"
fi

SORTMERNA_DIR="${DB_DIR}/sortmerna"
CHECKSUMS="${SORTMERNA_DIR}/SHA256SUMS"

log "database directory : ${DB_DIR}"
log "index directory    : ${SORTMERNA_DIR}"

# --------------------------------------------------------- make checksums ----
if [ "${MAKE_CHECKSUMS}" = "true" ]; then
    [ -d "${SORTMERNA_DIR}" ] || fail "index directory not found: ${SORTMERNA_DIR}"
    present=()
    for f in "${INDEX_FILES[@]}"; do
        [ -f "${SORTMERNA_DIR}/${f}" ] && present+=("$f")
    done
    [ "${#present[@]}" -gt 0 ] || fail "no index files present in ${SORTMERNA_DIR}"
    log "writing checksums for ${#present[@]} file(s)"
    if [ "${DRY_RUN}" = "true" ]; then
        log "DRY-RUN: would write ${CHECKSUMS}"
        exit 0
    fi
    ( cd "${SORTMERNA_DIR}" && sha256sum "${present[@]}" > SHA256SUMS )
    log "written: ${CHECKSUMS}"
    printf '%s\n' "Upload these files as release assets:" "${present[@]/#/  }"
    exit 0
fi

# ------------------------------------------------------------ verify mode ----
verify_installed() {
    local bad=0 f
    for f in "${INDEX_FILES[@]}"; do
        if [ ! -s "${SORTMERNA_DIR}/${f}" ]; then
            printf '  MISSING  %s\n' "$f"
            bad=$((bad + 1))
        else
            printf '  present  %-24s %10s bytes\n' "$f" \
                "$(stat -L -c %s "${SORTMERNA_DIR}/${f}" 2>/dev/null || echo '?')"
        fi
    done
    return "$bad"
}

if [ "${CHECK_ONLY}" = "true" ]; then
    [ -d "${SORTMERNA_DIR}" ] || fail "index directory not found: ${SORTMERNA_DIR}"
    if verify_installed; then
        log "index complete"
        exit 0
    fi
    fail "the index is incomplete; rerun without --check to download it"
fi

# ------------------------------------------------------- already installed ---
if [ "${FORCE}" != "true" ] && [ -d "${SORTMERNA_DIR}" ]; then
    if verify_installed >/dev/null 2>&1; then
        log "index already installed and complete; pass --force to re-download"
        exit 0
    fi
fi

# ---------------------------------------------------------------- download ---
[ -n "${BASE_URL}" ] || fail "--url is required to download the index"
BASE_URL="${BASE_URL%/}"

# Pick a downloader. Fail early rather than half way through a 1.2 GB transfer.
if command -v curl >/dev/null 2>&1; then
    DL=(curl -fL --retry 3 --retry-delay 5 -o)
elif command -v wget >/dev/null 2>&1; then
    DL=(wget --tries=3 --waitretry=5 -O)
else
    fail "neither curl nor wget is available"
fi

FILES=("${INDEX_FILES[@]}")
[ "${WITH_FASTA}" = "true" ] && FILES+=("rrna_merged.fasta")

if [ "${DRY_RUN}" = "true" ]; then
    log "DRY-RUN: would create ${SORTMERNA_DIR} and download:"
    for f in "${FILES[@]}"; do
        printf '  %s/%s\n' "${BASE_URL}" "$f"
    done
    log "DRY-RUN: would verify against ${CHECKSUMS} and rename .partial files"
    exit 0
fi

mkdir -p "${SORTMERNA_DIR}"

# Checksums are optional: a release without SHA256SUMS still installs, but says so.
if [ -f "${CHECKSUMS}" ]; then
    log "checksums: ${CHECKSUMS}"
else
    log "checksums: not present; integrity will not be verified"
fi

failed=0
for f in "${FILES[@]}"; do
    target="${SORTMERNA_DIR}/${f}"
    if [ -s "${target}" ] && [ "${FORCE}" != "true" ]; then
        log "skip ${f} (already present)"
    else
        partial="${target}.partial"
        log "downloading ${f}"
        if ! "${DL[@]}" "${partial}" "${BASE_URL}/${f}"; then
            log "download failed: ${f}"
            failed=$((failed + 1))
            continue
        fi
        mv -f "${partial}" "${target}"
    fi
done

[ "${failed}" -eq 0 ] || fail "${failed} file(s) failed to download; rerun to retry"

# ------------------------------------------------------------------- verify ---
if [ -f "${CHECKSUMS}" ]; then
    log "verifying SHA256 checksums"
    if ! ( cd "${SORTMERNA_DIR}" && sha256sum -c SHA256SUMS --ignore-missing ); then
        fail "checksum verification failed; delete the bad file(s) and rerun"
    fi
fi

if ! verify_installed; then
    fail "the index is still incomplete after downloading"
fi

log "done: ${SORTMERNA_DIR}"
log "1-qc.sh reads this through BOWTIE2_INDEX=<db>/sortmerna/rrna_index"
