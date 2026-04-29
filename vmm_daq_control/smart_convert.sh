#!/bin/bash
# smart_convert.sh
# Watch a directory for new pcapng files and run: convertFile -f <file>
# Only converts files smaller than 1 GB
# Usage: ./smart_convert.sh {start|stop|status} <watch_dir> <interface>

set -uo pipefail

### === Configuration ===

WATCH_DIR=$2
INTERFACE=$3
CONVERT_CMD="convertFile -geo config/geometry_${INTERFACE}.json -bc 44.444 -tac 60 -save [[0],[],[]] -df SRS -f"
MAX_SIZE_BYTES=$((1 * 1024 * 1024 * 1024))      # 1 GB in bytes

### === Functions ===

do_action() {
    local file="$1"
    
    # Check if file exists
    if [ ! -f "$file" ]; then
        echo "ACTION: file not found: $file"
        return 1
    fi
    
    # Get file size in bytes
    local filesize=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    
    if [ -z "$filesize" ]; then
        echo "ACTION: could not determine size of: $file"
        return 1
    fi
    
    # Convert bytes to human-readable format for logging
    local size_mb=$((filesize / 1024 / 1024))
    local size_gb=$(echo "scale=2; $filesize / 1024 / 1024 / 1024" | bc)
    
    # Check if file is too large
    if [ "$filesize" -gt "$MAX_SIZE_BYTES" ]; then
        echo "ACTION: SKIPPING large file (${size_gb} GB > 1 GB): $file"
        return 0  # Not an error, just skipped
    fi
    
    echo "ACTION: file size OK (${size_mb} MB), running: $CONVERT_CMD $file"

    if eval "$CONVERT_CMD $file" ; then
        echo "ACTION: convert succeeded for: $file"
        return 0
    else
        echo "ACTION: convert FAILED for: $file"
        return 1
    fi
}

process_path() {
    local path="$1"
    # Only process .pcapng files (adjust pattern if needed)
}

### === Main ===

while IFS='|' read -r ev filepath; do
    echo "Event: $ev -> $filepath"
    case "$filepath" in
        *.pcapng) do_action "$filepath" ;;
        *) echo "Skipping non-pcapng: $filepath" ;;
    esac
done < <(inotifywait -m -e close_write --include "${INTERFACE}" --format '%e|%w%f' --quiet "$WATCH_DIR")
