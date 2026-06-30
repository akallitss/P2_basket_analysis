#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hv_live_monitor.py

Stand-alone reproduction of the CAEN HV live monitoring that runs inside the
DAQ (Cosmic_Bench_DAQ_Control). During a run, hv_control.py polls every HV
channel and appends one row per readout to <sub_run>/hv_monitor.csv:

    timestamp, <slot>:<ch> power, <slot>:<ch> v0, <slot>:<ch> vmon,
               <slot>:<ch> imon, ...   (repeated for every channel)

      power : channel on (1) / off (0)
      v0    : set voltage          [V]
      vmon  : measured voltage     [V]
      imon  : measured current     [uA]

The DAQ flask dashboard simply re-reads the tail of this CSV every few seconds
and draws two time-series charts (voltage and current, one trace per channel).
This script does the same thing without needing the DAQ stack:

  * --live : tail the growing CSV and refresh two panels (V & I vs time) every
             few seconds, exactly like the online monitor. Works whether it
             points at a local copy or the DAQ-mounted file as it grows.
  * static (default): plot the whole saved hv_monitor.csv for a finished run
             and write a PNG.

Examples
  # Replay / inspect a finished run (writes hv_monitor.png next to the CSV)
  python3 hv_live_monitor.py --csv /path/to/.../hv_monitor.csv

  # Live monitor a run that is currently being written
  python3 hv_live_monitor.py --csv /path/to/.../hv_monitor.csv --live

  # Live monitor, last 500 readouts, refresh every 2 s
  python3 hv_live_monitor.py --csv .../hv_monitor.csv --live --tail 500 --interval 2
"""

import os
import argparse

import matplotlib
import pandas as pd


# --------------------------------------------------------------------------- #
# CSV parsing
# --------------------------------------------------------------------------- #
def parse_hv_csv(csv_path, tail=None):
    """
    Read hv_monitor.csv and split it into per-channel measurements.

    Returns (times, channels) where
      times    : pandas Series of datetime timestamps
      channels : dict keyed by 'slot:ch' -> dict with 'power', 'v0', 'vmon',
                 'imon' Series (aligned with `times`).
    """
    df = pd.read_csv(csv_path)
    if tail is not None and tail > 0:
        df = df.tail(tail)

    # ISO8601 parses both the old 1 s ('...HH:MM:SS') and new sub-second
    # ('...HH:MM:SS.fff') hv_monitor.csv timestamp formats.
    times = pd.to_datetime(df['timestamp'], format='ISO8601', errors='coerce')

    channels = {}
    for col in df.columns:
        if not col.endswith(' vmon'):
            continue
        key = col[:-len(' vmon')]  # 'slot:ch'
        channels[key] = {
            'power': df.get(f'{key} power'),
            'v0': df.get(f'{key} v0'),
            'vmon': df[f'{key} vmon'],
            'imon': df.get(f'{key} imon'),
        }
    return times, channels


def channel_sort_key(key):
    """Sort 'slot:ch' numerically by (slot, channel)."""
    try:
        slot, ch = key.split(':')
        return (int(slot), int(ch))
    except ValueError:
        return (9999, 9999)


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def draw(ax_v, ax_i, times, channels):
    """(Re)draw the voltage and current panels from scratch."""
    import matplotlib.dates as mdates

    ax_v.clear()
    ax_i.clear()

    keys = sorted(channels, key=channel_sort_key)
    for key in keys:
        ch = channels[key]
        powered = ''
        if ch['power'] is not None and len(ch['power']):
            powered = '' if bool(ch['power'].iloc[-1]) else ' (off)'
        ax_v.plot(times, ch['vmon'], lw=1.2, label=f'{key}{powered}')
        if ch['imon'] is not None:
            ax_i.plot(times, ch['imon'], lw=1.2, label=f'{key}{powered}')

    ax_v.set_ylabel('Voltage  [V]')
    ax_v.set_title('HV monitor — measured voltage (vmon) vs time')
    ax_v.grid(True, alpha=0.3)
    ax_v.legend(loc='center left', bbox_to_anchor=(1.01, 0.5),
                fontsize=8, ncol=1, title='slot:ch')

    ax_i.set_ylabel('Current  [µA]')
    ax_i.set_xlabel('Time')
    ax_i.set_title('Measured current (imon) vs time')
    ax_i.grid(True, alpha=0.3)
    ax_i.legend(loc='center left', bbox_to_anchor=(1.01, 0.5),
                fontsize=8, ncol=1, title='slot:ch')

    for ax in (ax_v, ax_i):
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    for label in ax_i.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')


def print_status(times, channels):
    """Print the latest readout, mirroring hv_control.py's screen output."""
    if times is None or len(times) == 0:
        print('No data yet.')
        return
    last_t = times.iloc[-1]
    print(f"\nHV monitor  {last_t}")
    for key in sorted(channels, key=channel_sort_key):
        ch = channels[key]
        slot, c = key.split(':')
        power = bool(ch['power'].iloc[-1]) if ch['power'] is not None else None
        v0 = ch['v0'].iloc[-1] if ch['v0'] is not None else float('nan')
        vmon = ch['vmon'].iloc[-1]
        imon = ch['imon'].iloc[-1] if ch['imon'] is not None else float('nan')
        print(f"  Slot {slot:>2} Ch {c:>2}: "
              f"power={'on ' if power else 'off'}, "
              f"v set={v0:>7.2f}, v mon={vmon:>8.2f} V, i mon={imon:>9.4f} uA")


# --------------------------------------------------------------------------- #
# Modes
# --------------------------------------------------------------------------- #
def run_static(csv_path, out_path, tail, show):
    import matplotlib.pyplot as plt

    times, channels = parse_hv_csv(csv_path, tail=tail)
    fig, (ax_v, ax_i) = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    draw(ax_v, ax_i, times, channels)
    fig.suptitle(os.path.basename(os.path.dirname(csv_path)) or 'hv_monitor',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    print_status(times, channels)
    if out_path:
        fig.savefig(out_path, dpi=130, bbox_inches='tight')
        print(f'\nSaved {out_path}')
    if show:
        plt.show()
    else:
        plt.close(fig)


def run_live(csv_path, interval, tail):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    fig, (ax_v, ax_i) = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    fig.suptitle(f'Live HV monitor — {csv_path}', fontsize=10)

    def update(_frame):
        try:
            times, channels = parse_hv_csv(csv_path, tail=tail)
        except (FileNotFoundError, pd.errors.EmptyDataError):
            return
        if times is None or len(times) == 0:
            return
        draw(ax_v, ax_i, times, channels)
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        print_status(times, channels)

    update(0)
    # Keep a reference so the animation is not garbage-collected.
    fig._hv_anim = FuncAnimation(fig, update, interval=interval * 1000,
                                 cache_frame_data=False)
    print(f'Live monitoring {csv_path} (refresh every {interval}s). '
          f'Close the window to stop.')
    plt.show()


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description='Stand-alone CAEN HV live monitor from hv_monitor.csv.')
    ap.add_argument('--csv', required=True, help='Path to hv_monitor.csv.')
    ap.add_argument('--live', action='store_true',
                    help='Tail the CSV and refresh continuously (online monitor).')
    ap.add_argument('--interval', type=float, default=5.0,
                    help='Live refresh interval in seconds (default 5).')
    ap.add_argument('--tail', type=int, default=None,
                    help='Only use the most recent N readouts '
                         '(default: all in static mode, 1000 in live mode).')
    ap.add_argument('--out', default=None,
                    help='Static mode: output PNG path '
                         '(default: hv_monitor.png next to the CSV).')
    ap.add_argument('--no-show', action='store_true',
                    help='Static mode: do not open a window, only save the PNG.')
    args = ap.parse_args()

    if not os.path.isfile(args.csv):
        ap.error(f'CSV not found: {args.csv}')

    if args.live:
        # Pick the first interactive backend that imports on this machine.
        for backend in ('TkAgg', 'QtAgg', 'GTK3Agg', 'MacOSX'):
            try:
                matplotlib.use(backend)
                # 'as _plt' so we don't bind the name `matplotlib` locally in
                # main() (which would shadow the module-level import).
                import matplotlib.pyplot as _plt  # noqa: F401  (force backend)
                break
            except Exception:
                continue
        else:
            ap.error('No interactive matplotlib backend available for --live; '
                     'run without --live to save a PNG instead.')
        tail = args.tail if args.tail is not None else 1000
        run_live(args.csv, args.interval, tail)
    else:
        if args.no_show:
            matplotlib.use('Agg')
        out = args.out
        if out is None:
            out = os.path.join(os.path.dirname(os.path.abspath(args.csv)),
                               'hv_monitor.png')
        run_static(args.csv, out, args.tail, show=not args.no_show)


if __name__ == '__main__':
    main()
