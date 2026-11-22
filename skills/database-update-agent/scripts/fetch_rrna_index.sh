#!/usr/bin/env bash
# Fetch the bowtie2 rRNA index used by 1-qc.sh for rRNA depletion.
#
# The index is ~1.2 GB and cannot live in the git repository, so it is published
# as a release asset. Two layouts are supported and the URL decides which runs:
#
#   archive mode   the URL ends in .tar.gz or .tgz. One download, then the six
#                  .bt2 files are extracted into <db>/sortmerna/.
#   files mode     any other URL. Each of the six .bt2 files is fetched
#                  separately, which lets you replace one corrupt file without
#                  pulling 1.2 GB again.
#
# Only the .bt2 files are needed. bowtie2 is invoked with -x <prefix>, which
# reads the index and never touches the source FASTA. rrna_merged.fasta is
# therefore published optionally; pass --with-fasta if you also want it for
# rebuilding the index later. That flag applies to files mode only.
#
# Usage
#   bash fetch_rrna_index.sh --db-dir <DB>            # default release asset
#   bash fetch_rrna_index.sh --url <url>              # explicit source
#   bash fetch_rrna_index.sh --db-dir <DB> --check    # verify only
#   bash fetch_rrna_index.sh --db-dir <DB> --make-checksums
#
# Exit codes
#   0  success, or the index is already installed and intact
#   1  a download, checksum, or installation step failed
#   2  bad arguments or an unresolved database directory
#
# Every failure is reported on stderr with the path involved. Nothing you might
# need is deleted: a partial download stays as .partial, and a failed extraction
# stays under .extract-<random> so it can be inspected.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_NAMES="config.yaml config-cell-paper.yaml cell_submit-paper.yaml"

# Default release asset, used when --url is omitted so a fresh clone needs only
# --db-dir.
DEFAULT_URL="https://github.com/Yu-Lab-Genomics/Hanta-Agent/releases/download/sortmerna-v1.0.0/rrna_index.tar.gz"

# SHA256 of that archive, used in archive mode only. Leave empty to skip the
# check and print a warning. After uploading a new archive, run
# --make-checksums to obtain the value and paste it here. A value in config.yaml
# (bowtie2.fetch_sha256) or on the command line (--sha256) takes precedence.
BUILTIN_ARCHIVE_SHA256=""

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
KEEP_ARCHIVE="false"
ARCHIVE_SHA256=""
BASE_URL=""
DB_DIR=""
CONFIG_YAML=""

usage() {
    cat <<'EOF'
Usage:
  bash fetch_rrna_index.sh --db-dir <path> [options]

Options:
  --url <url>          Release asset to fetch. Defaults to the bundled archive:
                       .../releases/download/sortmerna-v1.0.0/rrna_index.tar.gz
                       A URL ending in .tar.gz or .tgz selects archive mode; any
                       other URL selects files mode and is used as a base URL.
  --db-dir <path>      Database directory. Defaults to <project_dir>/DB, read
                       from a configuration file found upward from this script.
  -c, --config <path>  Explicit configuration file.
  --sha256 <hex>       Expected SHA256 of the archive. Archive mode only;
                       overrides the value built into the script.
  --keep-archive       Keep the .tar.gz after extraction. Archive mode only;
                       by default it is deleted to save about 1 GB.
  --check              Verify an already installed index and exit. Writes nothing.
  --make-checksums     Write SHA256SUMS for the files present in <db>/sortmerna,
                       and print the SHA256 of any archive found there. Run this
                       before uploading the assets.
  --with-fasta         Files mode only. Also download rrna_merged.fasta (needed
                       only to rebuild the index; the pipeline never reads it).
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
        --keep-archive)   KEEP_ARCHIVE="true"; shift ;;
        --sha256)         ARCHIVE_SHA256="${2:-}"; shift 2 ;;
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

# bowtie2.fetch_url and bowtie2.fetch_sha256 in config.yaml apply here too, so the
# index location is defined in one place. Precedence for both values is:
#   command line  >  config.yaml  >  the constants at the top of this script.
read_bowtie2_key() {
    awk -v k="$1" '
        /^[[:space:]]*#/ { next }
        /^bowtie2[[:space:]]*:/ { inb = 1; next }
        inb && /^[^[:space:]]/ { inb = 0 }
        inb {
            line = $0
            sub(/^[[:space:]]*/, "", line)
            if (line ~ "^" k "[[:space:]]*:") {
                sub("^" k "[[:space:]]*:[[:space:]]*", "", line)
                sub(/[[:space:]]+#.*$/, "", line)
                gsub(/"/, "", line)
                gsub(/\047/, "", line)
                print line
                exit
            }
        }' "$2" 2>/dev/null
}

if [ -f "${CONFIG_YAML}" ]; then
    _cfg_url="$(read_bowtie2_key fetch_url "${CONFIG_YAML}")"
    _cfg_sha="$(read_bowtie2_key fetch_sha256 "${CONFIG_YAML}")"
    [ -n "${BASE_URL}" ] || [ -z "${_cfg_url}" ] || BASE_URL="${_cfg_url}"
    [ -n "${ARCHIVE_SHA256}" ] || [ -z "${_cfg_sha}" ] || ARCHIVE_SHA256="${_cfg_sha}"
fi
[ -n "${ARCHIVE_SHA256}" ] || ARCHIVE_SHA256="${BUILTIN_ARCHIVE_SHA256}"

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
    # An archive in the same directory is what archive mode downloads, so its
    # hash is what ARCHIVE_SHA256 or --sha256 needs.
    for a in "${SORTMERNA_DIR}"/*.tar.gz "${SORTMERNA_DIR}"/*.tgz; do
        [ -f "$a" ] || continue
        log "archive SHA256; paste into ARCHIVE_SHA256 or pass via --sha256:"
        sha256sum "$a"
    done
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
# An explicit --url wins, then config.yaml's bowtie2.fetch_url, then the bundled
# default. The first two were already resolved just above.
[ -n "${BASE_URL}" ] || BASE_URL="${DEFAULT_URL}"

# Pick a downloader. Fail early rather than half way through a 1.2 GB transfer.
if command -v curl >/dev/null 2>&1; then
    DL=(curl -fL --retry 3 --retry-delay 5 -o)
elif command -v wget >/dev/null 2>&1; then
    DL=(wget --tries=3 --waitretry=5 -O)
else
    fail "neither curl nor wget is available"
fi

is_archive_url() {
    case "$1" in
        *.tar.gz|*.tgz) return 0 ;;
        *)              return 1 ;;
    esac
}

# ------------------------------------------------------------- archive mode ---
# One downloaded .tar.gz holding all six index files. The archive may store them
# at its root, which is the usual case, or inside a directory; either way each
# file is located by name and moved up, so no assumption is made about the top
# level and a stray wrapper directory is harmless.
install_from_archive() {
    local url="$1"
    local name archive partial tmp f src
    name="$(basename "${url}")"
    archive="${SORTMERNA_DIR}/${name}"

    if [ "${DRY_RUN}" = "true" ]; then
        log "DRY-RUN: archive mode"
        log "DRY-RUN: source  ${url}"
        log "DRY-RUN: archive ${archive}"
        log "DRY-RUN: target  ${SORTMERNA_DIR}/ <six .bt2 files>"
        if [ "${KEEP_ARCHIVE}" = "true" ]; then
            log "DRY-RUN: archive would be kept"
        else
            log "DRY-RUN: archive would be deleted after extraction"
        fi
        return 0
    fi

    mkdir -p "${SORTMERNA_DIR}"

    if [ -s "${archive}" ] && [ "${FORCE}" != "true" ]; then
        log "reusing archive already present: ${name}"
    else
        partial="${archive}.partial"
        log "downloading ${url}"
        if ! "${DL[@]}" "${partial}" "${url}"; then
            fail "download failed: ${url} (partial data left at ${partial})"
        fi
        mv -f "${partial}" "${archive}"
    fi

    if [ -n "${ARCHIVE_SHA256}" ]; then
        log "verifying archive SHA256"
        local actual
        actual="$(sha256sum "${archive}" | awk '{print $1}')"
        if [ "${actual}" != "${ARCHIVE_SHA256}" ]; then
            fail "archive checksum mismatch: expected ${ARCHIVE_SHA256}, got ${actual}. Delete ${archive} and rerun, or pass the correct value with --sha256"
        fi
        log "archive checksum OK"
    else
        log "no SHA256 configured for the archive; integrity not verified"
    fi

    # Extraction needs roughly as much space again as the archive itself. This is
    # a warning rather than a failure because the ratio depends on compression.
    local avail_bytes archive_bytes
    avail_bytes="$(df -Pk "${SORTMERNA_DIR}" | awk 'NR==2 {print $4 * 1024}')"
    archive_bytes="$(stat -L -c %s "${archive}" 2>/dev/null || echo 0)"
    if [ -n "${avail_bytes}" ] && [ "${archive_bytes}" -gt 0 ] \
       && [ "${avail_bytes}" -lt "$((archive_bytes * 2))" ]; then
        log "warning: free space under ${SORTMERNA_DIR} is below twice the archive size; extraction may fail"
    fi

    tmp="$(mktemp -d "${SORTMERNA_DIR}/.extract-XXXXXX")" \
        || fail "cannot create a temporary directory under ${SORTMERNA_DIR}"

    log "extracting ${name}"
    if ! tar -xzf "${archive}" -C "${tmp}"; then
        fail "extraction failed; extracted leftovers kept in ${tmp}"
    fi

    local -a missing=()
    for f in "${INDEX_FILES[@]}"; do
        src="$(find "${tmp}" -type f -name "${f}" -print -quit)"
        [ -n "${src}" ] || missing+=("${f}")
    done
    if [ "${#missing[@]}" -gt 0 ]; then
        fail "archive does not contain: ${missing[*]} (extracted tree kept in ${tmp})"
    fi

    for f in "${INDEX_FILES[@]}"; do
        src="$(find "${tmp}" -type f -name "${f}" -print -quit)"
        mv -f "${src}" "${SORTMERNA_DIR}/${f}"
    done
    rm -rf "${tmp}"
    log "installed ${#INDEX_FILES[@]} index files into ${SORTMERNA_DIR}"

    if [ "${KEEP_ARCHIVE}" = "true" ]; then
        log "archive kept: ${archive}"
    else
        rm -f "${archive}"
        log "archive removed: ${name} (pass --keep-archive to keep it)"
    fi
}

# --------------------------------------------------------------- files mode ---
# Six separate downloads against a base URL. Kept because it lets you replace a
# single corrupt file without pulling the whole 1.2 GB again.
install_from_files() {
    local url="$1"
    local -a files=("${INDEX_FILES[@]}")
    [ "${WITH_FASTA}" = "true" ] && files+=("rrna_merged.fasta")

    if [ "${DRY_RUN}" = "true" ]; then
        log "DRY-RUN: files mode"
        log "DRY-RUN: would create ${SORTMERNA_DIR} and download:"
        for f in "${files[@]}"; do
            printf '  %s/%s\n' "${url}" "$f"
        done
        log "DRY-RUN: would verify against ${CHECKSUMS} and rename .partial files"
        return 0
    fi

    mkdir -p "${SORTMERNA_DIR}"

    # Checksums are optional: a release without SHA256SUMS still installs, but says so.
    if [ -f "${CHECKSUMS}" ]; then
        log "checksums: ${CHECKSUMS}"
    else
        log "checksums: not present; integrity will not be verified"
    fi

    local failed=0 f target partial
    for f in "${files[@]}"; do
        target="${SORTMERNA_DIR}/${f}"
        if [ -s "${target}" ] && [ "${FORCE}" != "true" ]; then
            log "skip ${f} (already present)"
        else
            partial="${target}.partial"
            log "downloading ${f}"
            if ! "${DL[@]}" "${partial}" "${url}/${f}"; then
                log "download failed: ${f}"
                failed=$((failed + 1))
                continue
            fi
            mv -f "${partial}" "${target}"
        fi
    done

    [ "${failed}" -eq 0 ] || fail "${failed} file(s) failed to download; rerun to retry"

    if [ -f "${CHECKSUMS}" ]; then
        log "verifying SHA256 checksums"
        if ! ( cd "${SORTMERNA_DIR}" && sha256sum -c SHA256SUMS --ignore-missing ); then
            fail "checksum verification failed; delete the bad file(s) and rerun"
        fi
    fi
}

# -------------------------------------------------------------------- drive ---
if is_archive_url "${BASE_URL}"; then
    if [ "${WITH_FASTA}" = "true" ]; then
        log "warning: --with-fasta has no effect in archive mode; the archive defines its own contents"
    fi
    install_from_archive "${BASE_URL}"
else
    BASE_URL="${BASE_URL%/}"
    install_from_files "${BASE_URL}"
fi

[ "${DRY_RUN}" = "true" ] && exit 0

# ---------------------------------------------------------------- verify ------
if ! verify_installed; then
    fail "the index is still incomplete after downloading"
fi

log "done: ${SORTMERNA_DIR}"
log "1-qc.sh reads this through BOWTIE2_INDEX=<db>/sortmerna/rrna_index"
