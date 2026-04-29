#!/bin/bash
# convert.sh
# Watch a directory for new pcapng files and run: convertFile -f <file>
# Usage: ./convert.sh {start|stop|status}

set -uo pipefail

### === Configuration ===

WATCH_DIR=$2                                     # directory to watch (change as needed)
INTERFACE=$3
CONVERT_CMD="convertFile -geo config/geometry_${INTERFACE}.json -bc 44.444 -tac 60 -save [[0],[],[]] -df SRS -f"

### === Functions ===

do_action() {
    local file="$1"
    echo "ACTION: running: $CONVERT_CMD $file"
    # run convert command; ensure proper quoting and preserve exit code
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

