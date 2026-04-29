#!/bin/bash
# smart_convert.sh
# Watch a directory for new pcapng files and run: convertFile -f <file>
# Only converts files smaller than 1 GB
# Usage: ./smart_convert.sh {start|stop|status} <watch_dir> <interface>

set -uo pipefail

### === Configuration ===

WATCH_DIR=$2
INTERFACE=$3
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONVERT_CMD="convertFile \
    -geo ${SCRIPT_DIR}/config/geometry_${INTERFACE}.json \
    -bc 44.444 -tac 60 \
    -save [[0],[],[]] \
    -df SRS -f"

MAX_SIZE_BYTES=$((1 * 1024 * 1024 * 1024))      # 1 GB in bytes
MEM_LIMIT="4G"                                 # kill if exceeds 4GB RAM
LOCK_FILE="/tmp/convert_lock_${INTERFACE}"

### === Functions ===

do_action() {
    local file="$1"

    # Check file exists
    if [ ! -f "$file" ]; then
        echo "[WARN] file not found: $file"
        return 1
    fi

    # Get file size
    local filesize
    filesize=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    if [ -z "$filesize" ]; then
        echo "[WARN] could not determine size of: $file"
        return 1
    fi

    local size_mb=$(( filesize / 1024 / 1024 ))
    local size_gb
    size_gb=$(echo "scale=2; $filesize / 1024 / 1024 / 1024" | bc)

    # Skip files that are too large
    if [ "$filesize" -gt "$MAX_SIZE_BYTES" ]; then
        echo "[SKIP] file too large (${size_gb} GB > 1 GB): $file"
        return 0
    fi

    # Prevent parallel conversions for same interface
    while [ -f "$LOCK_FILE" ]; do
        echo "[WAIT] Another conversion in progress, waiting... ($file)"
        sleep 5
    done

    # Acquire lock
    touch "$LOCK_FILE"
    trap 'rm -f "$LOCK_FILE"; trap - RETURN' RETURN

    echo "[INFO] Converting (${size_mb} MB): $file"

    # Run with memory cap — kills gracefully instead of crashing machine
    if systemd-run --scope \
            -p MemoryMax="$MEM_LIMIT" \
            -p MemorySwapMax=0 \
            -- bash -c "$CONVERT_CMD $file"; then
        echo "[OK] Converted: $file"
    else
        echo "[ERROR] Conversion failed (exit $?): $file"
        rm -f "$LOCK_FILE"
        return 1
    fi

    rm -f "$LOCK_FILE"
    sleep 2  # let RAM settle before next conversion
}

### === Main watch loop ===
echo "[START] Watching $WATCH_DIR for interface $INTERFACE..."

inotifywait \
    -m \
    -e close_write \
    --format '%e|%w%f' \
    --quiet \
    "$WATCH_DIR" \
| while IFS='|' read -r ev filepath; do
    echo "[EVENT] $ev -> $filepath"
    case "$filepath" in
        *${INTERFACE}*.pcapng)
            do_action "$filepath" ;;
        *.pcapng)
            echo "[SKIP] pcapng but wrong interface: $filepath" ;;
        *)
            echo "[SKIP] not a pcapng: $filepath" ;;
    esac
done
