#!/bin/bash

set -euo pipefail

### Configuration ###

#DURATION=600
DURATION=44
OUTBASE=/local/p2/p2data/cern_202511
WAIT_TIMEOUT=60
SLEEP_INTERVAL=5

LOGDIR=logs
mkdir -p "$LOGDIR"

usage() { echo "Usage: $0 {start|stop|status} <run_name>"; exit 1; }

ACTION=${1:-}
RUN=${2:-}
[ -z "$ACTION" ] && usage
[ -z "$RUN" ] && [ "$ACTION" == "start" ] && usage

RUNDIR=${OUTBASE}/${RUN:-"latest"}
mkdir -p "$RUNDIR"

PID_DAQ_SRS=/tmp/pid_daq_srs
PID_DAQ_ALINX=/tmp/pid_daq_alinx
PID_CONVERT_SRS=/tmp/pid_convert_srs
PID_CONVERT_ALINX=/tmp/pid_convert_alinx

### Functions ###

start() {

    echo "[START] SRS: Starting DAQ..."
    ./loop_daq.sh enx0c37968d3d99 "$DURATION" "$RUNDIR" &
    echo $! > "$PID_DAQ_SRS"

    echo "[START] ALINX: Checking status..."
    alinx-sc --config-file config/config_alinx_noThresholds.json --read-link-status
    echo "[START] ALINX: Switching acquisition ON..."
    alinx-sc --config-file config/config_alinx_noThresholds.json --acq-on

    echo "[START] ALINX: Starting DAQ..."
    ./loop_daq.sh enp4s0f1 "$DURATION" "$RUNDIR" &
    echo $! > "$PID_DAQ_ALINX"

    echo "[START] File conversion watchers..."
    ./smart_convert.sh start "$RUNDIR" enx0c37968d3d99 > "$LOGDIR/convert_enx0c37968d3d99.log" 2>&1 &
    echo $! > "$PID_CONVERT_SRS"

    ./smart_convert.sh start "$RUNDIR" enp4s0f1 > "$LOGDIR/convert_enp4s0f1.log" 2>&1 &
    echo $! > "$PID_CONVERT_ALINX"

    echo "[START] All started."

}

stop() {
    echo "[STOP] Stopping DAQs..."
    ./stop_process.sh "ALINX DAQ" "$PID_DAQ_ALINX"
    ./stop_process.sh "SRS DAQ"   "$PID_DAQ_SRS"

    # ── Wait for conversions only if watchers are actually running ──
    local srs_pid alinx_pid srs_running alinx_running
    srs_pid=$(cat "$PID_CONVERT_SRS" 2>/dev/null || echo "")
    alinx_pid=$(cat "$PID_CONVERT_ALINX" 2>/dev/null || echo "")

    srs_running=false
    alinx_running=false
    [ -n "$srs_pid" ]   && kill -0 "$srs_pid"   2>/dev/null && srs_running=true || true
    [ -n "$alinx_pid" ] && kill -0 "$alinx_pid" 2>/dev/null && alinx_running=true || true

    if $srs_running || $alinx_running; then
        echo "[INFO] Conversion watchers running, waiting for pending files..."
        local npcap nroot remaining
        local end_time=$(( SECONDS + WAIT_TIMEOUT ))

        while [ $SECONDS -lt $end_time ]; do
            npcap=$(find "$RUNDIR" -maxdepth 1 -type f -name "*.pcapng" | wc -l)
            nroot=$(find "$RUNDIR" -maxdepth 1 -type f -name "*.root"   | wc -l)
            remaining=$(( npcap - nroot ))

            if [ "$remaining" -le 0 ]; then
                echo "[OK] All pcapng files converted."
                sleep 5  # let last conversion fully flush to disk
                break
            fi
            echo "[INFO] $remaining file(s) still pending. Waiting..."
            sleep "$SLEEP_INTERVAL"
        done

        # Warn if timeout reached
        npcap=$(find "$RUNDIR" -maxdepth 1 -type f -name "*.pcapng" | wc -l)
        nroot=$(find "$RUNDIR" -maxdepth 1 -type f -name "*.root"   | wc -l)
        if [ $(( npcap - nroot )) -gt 0 ]; then
            echo "[WARN] Timeout: $(( npcap - nroot )) file(s) not yet converted."
            echo "[WARN] Convert manually later with:"
            echo "[WARN]   for f in $RUNDIR/*.pcapng; do convertFile ... -f \$f; done"
        fi
    else
        echo "[INFO] Conversion watchers not running — skipping wait."
        echo "[INFO] Unconverted pcapng files (if any) in: $RUNDIR"
    fi

    echo "[STOP] Stopping watchers..."
    ./stop_process.sh "SRS File Conversion Watcher"   "$PID_CONVERT_SRS"
    ./stop_process.sh "ALINX File Conversion Watcher" "$PID_CONVERT_ALINX"
    pkill -f inotifywait || true  # don't fail if already dead

    echo "[STOP] ALINX: Switching acquisition OFF..."
    alinx-sc --config-file config/config_alinx_noThresholds.json --acq-off
    alinx-sc --config-file config/config_alinx_noThresholds.json --read-link-status

    echo "[STOP] All processes stopped."
}

status() {

    for p in "$PID_DAQ_SRS" "$PID_DAQ_ALINX" "$PID_CONVERT_SRS" "$PID_CONVERT_ALINX"; do
        if [ -f "$p" ]; then
            pid=$(cat "$p")
            if kill -0 "$pid" 2>/dev/null; then
                echo "[RUNNING] $(basename "$p"): PID $pid"
            else
                echo "[STALE]   $(basename "$p"): PID $pid (not alive)"
            fi
        else
            echo "[STOPPED] $(basename "$p")"
        fi
    done

}

### Dispatch ###

case "$ACTION" in
    start)  start ;;
    stop)   stop ;;
    status) status ;;
    *)      usage ;;
esac
