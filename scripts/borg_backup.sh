#!/usr/bin/env bash
#
# BACK UP PROJECT FILES
#
# Only includes things like logs and R objects, not all intermediate big results
# files such as BAM, FASTQ and the likes. Requires an environment variable with
# the path to the borg repository as well as a pattern file for which files to
# ignore/include.
#
# Usage:
#   For dry run and list of files that would be backed up:
#       `bash borg_backup.sh`
#
#   For the actual backup:
#       `bash borg_backup.sh run`
#
#
# NOTE for rkn servers: borg is only usable via docker container/distrobox.
#
# Markus Riedl
# Date: 2024-09-11
#

if [[ "$1" == "-n" ]]; then
    dry="-n"
fi

BACKUP_DIRS=(results resources plots)

borg create \
    ${dry:-"--stat"} \
    ${BORG_BACKUP_REPO}::$(date -Iminutes) \
    "${BACKUP_DIRS[@]}"

# ================================================
#   Semi old script which was still too complicated
# ================================================

: <<'END_COMMENT'


usage() {
    echo "Usage: ${0} [-m <message>] <files...>"
}

[[ "$1" == "-h" ]] && usage && exit 1

if [[ "$1" == "-m" ]]; then
    shift
    message="${1}"
    shift
elif [[ "$1" == "-n" ]]; then
    shift
    dry="-n"
fi

#
# Backup
#

borg create \
    --stat --list ${dry:-} \
    -e work -e .nextflow -e archive -e data -e renv/library \
    ${BORG_BACKUP_REPO}::$(date -Iminutes)_$(git rev-parse --short HEAD)_${message:-backup} \
    .

END_COMMENT

# ================================================
#
#
# OLD SCRIPT THAT BACKED UP EVERYTHING
#
: <<'END_COMMENT'
set -euo pipefail

if [[ "$#" -gt 0 && "$1" == "run" ]]; then
    args="--stats --list --filter AME"
else
    args="--dry-run --list --filter -"
fi

borg \
    create \
    ${args} \
    --patterns-from .borgignore \
    ${BORG_BACKUP_REPO}::$(date -Iminutes)_$(git rev-parse --short HEAD) \
    .
END_COMMENT
