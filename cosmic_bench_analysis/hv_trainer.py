#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hv_trainer.py

Software HV "training" / conditioning controller that replaces the CAEN
crate's hard over-current trip with a soft back-off-and-re-ramp loop, so a
sparky detector can be conditioned instead of having its HV killed.

Background
----------
A CAEN channel goes into OVC (over-current) when imon reaches its current
limit ISet; in that state it runs constant-current and vmon sags. If OVC
persists longer than the channel TRIP time (e.g. 30 s) the crate kills the
channel. During conditioning a fresh/sparky detector spends a lot of time in
OVC, so the crate keeps killing the HV. We want it to instead back the voltage
off for a moment, let the current settle, and slowly climb back toward the
target -- the manual "training" procedure, automated.

Controller (state machine)
---------------------------
Per readout (timestamp, vmon, imon) the controller updates a commanded set
voltage `vset` that would be pushed to the crate (caen.set_ch_v0):

  RAMP_UP : current safe (imon < i_safe) and held safe for `dwell` s
            -> step vset up by `v_step_up` toward `v_target`.
  BACKOFF : OVC (imon >= i_comp_frac*ICOMP) held for >= `backoff_after` s
            -> drop vset by `v_step_down` (>= v_floor), reset timers, and
               enter a recovery cooldown of `recover_dwell` s before any
               further ramp-up. backoff_after is set well below the crate
               TRIP time, so the crate never reaches its kill condition.
  HOLD    : otherwise keep vset.

`backoff_after` << TRIP time is the whole point: the longest continuous OVC
the trainer ever tolerates is ~backoff_after, so the crate's 30 s kill never
fires.

This file provides:
  * TrainingController  -- pure decision logic, no hardware dependency.
  * backtest_csv()      -- replay a recorded hv_monitor.csv channel through the
                           controller and plot/summarise it against the real
                           trace and the crate's would-be kill timeline.
  * run_live()          -- drive a real CAEN channel (needs caen_hv_py, guarded
                           import); same set_ch_v0 calls as hv_control.py.

IMPORTANT -- the backtest is OPEN LOOP. We have the recorded imon(t) but cannot
know how the detector current would have responded to the voltages the
controller commands. The backtest therefore validates the *decision logic and
timing* against real spark dynamics (when would the crate kill? when would the
trainer intervene, and by how much margin?). It is not a closed-loop physics
simulation of the conditioning.

Usage
  # Backtest the sparky channel from the overnight run
  python3 hv_trainer.py --csv /path/to/hv_monitor.csv --channel 1:0

  # Scan a more aggressive back-off and a lower target
  python3 hv_trainer.py --csv .../hv_monitor.csv --channel 1:0 \
      --backoff-after 4 --v-target 420 --v-step-down 30
"""

import os
import argparse
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# Controller
# --------------------------------------------------------------------------- #
@dataclass
class TrainerConfig:
    i_comp: float = 10.0        # crate current limit ISet [uA] -> OVC at this
    i_comp_frac: float = 0.95   # treat imon >= frac*i_comp as OVC
    i_safe: float = 2.0         # imon below this is "settled" [uA]
    v_target: float = 420.0     # conditioning target [V]
    v_floor: float = 0.0        # never command below this [V]
    v_step_up: float = 10.0     # ramp-up step [V]
    v_step_down: float = 25.0   # back-off drop per OVC event [V]
    dwell: float = 30.0         # safe-current hold before stepping up [s]
    recover_dwell: float = 30.0  # cooldown after a back-off before ramping [s]
    backoff_after: float = 6.0  # continuous OVC before backing off [s]
    caen_trip_time: float = 30.0  # crate's OVC->kill time, for comparison [s]


@dataclass
class TrainingController:
    cfg: TrainerConfig
    vset: float = None                 # current commanded set voltage [V]
    _ovc_run: float = 0.0              # continuous OVC duration [s]
    _safe_run: float = 0.0             # continuous settled duration [s]
    _cooldown: float = 0.0             # remaining recovery cooldown [s]
    events: list = field(default_factory=list)  # (t, kind, vset, imon)

    def __post_init__(self):
        if self.vset is None:
            self.vset = self.cfg.v_target

    def step(self, t, imon, dt):
        """Advance one readout. Returns (vset, state)."""
        c = self.cfg
        ovc = imon >= c.i_comp_frac * c.i_comp
        safe = imon < c.i_safe

        self._ovc_run = self._ovc_run + dt if ovc else 0.0
        self._safe_run = self._safe_run + dt if safe else 0.0
        self._cooldown = max(0.0, self._cooldown - dt)

        state = 'HOLD'
        if self._ovc_run >= c.backoff_after:
            # Soft trip: back off before the crate would kill.
            new = max(self.vset - c.v_step_down, c.v_floor)
            self.events.append((t, 'BACKOFF', new, imon))
            self.vset = new
            self._ovc_run = 0.0
            self._safe_run = 0.0
            self._cooldown = c.recover_dwell
            state = 'BACKOFF'
        elif (self.vset < c.v_target and self._safe_run >= c.dwell
              and self._cooldown <= 0.0):
            new = min(self.vset + c.v_step_up, c.v_target)
            if new != self.vset:
                self.events.append((t, 'STEP_UP', new, imon))
            self.vset = new
            self._safe_run = 0.0
            state = 'RAMP_UP'
        return self.vset, state


# --------------------------------------------------------------------------- #
# Crate kill timeline (what the CAEN 30 s rule would do on this trace)
# --------------------------------------------------------------------------- #
def crate_kill_events(secs, imon, power, cfg):
    """Return list of (t_idx) where continuous OVC first reaches trip time."""
    thr = cfg.i_comp_frac * cfg.i_comp
    kills, run, start = [], 0.0, None
    for k in range(len(secs)):
        ovc = (imon[k] >= thr) and (power[k] == 1)
        if ovc:
            if start is None:
                start = k
            run = secs[k] - secs[start]
            if run >= cfg.caen_trip_time:
                kills.append(k)
                run, start = 0.0, None  # crate would have to be re-enabled
        else:
            run, start = 0.0, None
    return kills


# --------------------------------------------------------------------------- #
# Backtest
# --------------------------------------------------------------------------- #
def backtest_csv(csv_path, channel, cfg, out_path, show, v_start=None):
    import matplotlib
    if not show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    df = pd.read_csv(csv_path)
    # ISO8601 parses both the old 1 s and new sub-second timestamp formats.
    t = pd.to_datetime(df['timestamp'], format='ISO8601')
    vmon = df[f'{channel} vmon'].values
    imon = df[f'{channel} imon'].values
    power = df[f'{channel} power'].values
    secs = (t - t.iloc[0]).dt.total_seconds().values

    # Run controller over the recorded current trace.
    ctrl = TrainingController(cfg, vset=(v_start if v_start is not None
                                         else cfg.v_target))
    vset_cmd, states = np.empty(len(df)), []
    for k in range(len(df)):
        dt = secs[k] - secs[k - 1] if k > 0 else 0.0
        vs, st = ctrl.step(t.iloc[k], imon[k], dt)
        vset_cmd[k] = vs
        states.append(st)

    kills = crate_kill_events(secs, imon, power, cfg)
    backoffs = [e for e in ctrl.events if e[1] == 'BACKOFF']
    stepups = [e for e in ctrl.events if e[1] == 'STEP_UP']

    # ---- console summary ----
    thr = cfg.i_comp_frac * cfg.i_comp
    ovc_mask = (imon >= thr) & (power == 1)
    print(f"\n=== Backtest {channel}  (open-loop replay of recorded imon) ===")
    print(f"  compliance ISet={cfg.i_comp} uA, OVC thr={thr:.1f} uA, "
          f"crate TRIP={cfg.caen_trip_time:.0f}s, trainer back-off@{cfg.backoff_after:.0f}s")
    print(f"  OVC samples: {int(ovc_mask.sum())} / {int((power==1).sum())} powered")
    if kills:
        print(f"  CRATE WOULD KILL HV {len(kills)}x; first at "
              f"{t.iloc[kills[0]]} (vmon~{vmon[kills[0]]:.0f} V)")
    else:
        print("  crate would NOT reach its 30 s kill on this trace")
    print(f"  trainer back-offs: {len(backoffs)} ; ramp-up steps: {len(stepups)}")
    if backoffs:
        print(f"  first trainer back-off at {backoffs[0][0]} "
              f"-> commanded {backoffs[0][2]:.0f} V "
              f"({'before' if not kills or backoffs[0][0] < t.iloc[kills[0]] else 'after'} crate kill)")
    longest_ovc = 0.0
    run = 0.0
    for k in range(len(secs)):
        if ovc_mask[k]:
            run += (secs[k] - secs[k-1]) if k > 0 else 0.0
            longest_ovc = max(longest_ovc, run)
        else:
            run = 0.0
    print(f"  longest continuous OVC in data: {longest_ovc:.0f}s "
          f"(trainer caps interventions at ~{cfg.backoff_after:.0f}s)")

    # ---- plot ----
    fig, (axv, axi) = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

    # shade OVC regions
    in_ovc = False
    for k in range(len(secs)):
        if ovc_mask[k] and not in_ovc:
            x0, in_ovc = t.iloc[k], True
        elif not ovc_mask[k] and in_ovc:
            for ax in (axv, axi):
                ax.axvspan(x0, t.iloc[k], color='red', alpha=0.08, lw=0)
            in_ovc = False

    axv.plot(t, vmon, lw=1.0, color='C0', label='vmon (recorded)')
    axv.plot(t, vset_cmd, lw=1.4, color='C2', ls='--',
             label='vset commanded by trainer')
    axv.axhline(cfg.v_target, color='gray', lw=0.8, ls=':', label='target')
    for e in backoffs:
        axv.axvline(e[0], color='darkorange', lw=0.8, alpha=0.6)
    for k in kills:
        axv.axvline(t.iloc[k], color='red', lw=1.5)
        axv.annotate('CAEN kill', (t.iloc[k], cfg.v_target),
                     color='red', fontsize=8, rotation=90, va='top')
    axv.set_ylabel('Voltage [V]')
    axv.set_title(f'HV trainer backtest — channel {channel} '
                  f'(open-loop replay; red bands = OVC, orange = trainer back-off)')
    axv.grid(True, alpha=0.3)
    axv.legend(loc='lower right', fontsize=8)

    axi.plot(t, imon, lw=1.0, color='C3', label='imon (recorded)')
    axi.axhline(thr, color='red', lw=0.8, ls=':', label=f'OVC thr {thr:.1f} µA')
    axi.axhline(cfg.i_safe, color='green', lw=0.8, ls=':',
                label=f'i_safe {cfg.i_safe:.1f} µA')
    for e in backoffs:
        axi.axvline(e[0], color='darkorange', lw=0.8, alpha=0.6)
    axi.set_ylabel('Current [µA]')
    axi.set_xlabel('Time')
    axi.grid(True, alpha=0.3)
    axi.legend(loc='upper right', fontsize=8)
    axi.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    for lab in axi.get_xticklabels():
        lab.set_rotation(30); lab.set_horizontalalignment('right')

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=130, bbox_inches='tight')
        print(f"\nSaved {out_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Live driver (real crate) — guarded import, mirrors hv_control.py
# --------------------------------------------------------------------------- #
def run_live(ip, username, password, slot, channel, cfg, poll=2.0):
    import time
    from caen_hv_py.CAENHVController import CAENHVController  # noqa

    print(f"Live training slot {slot} ch {channel} -> target {cfg.v_target} V")
    with CAENHVController(ip, username, password) as caen:
        ctrl = TrainingController(cfg, vset=caen.get_ch_v0(slot, channel))
        if not caen.get_ch_power(slot, channel):
            caen.set_ch_pw(slot, channel, 1)
        last = time.time()
        while ctrl.vset < cfg.v_target or True:  # run until interrupted
            imon = caen.get_ch_imon(slot, channel)
            vmon = caen.get_ch_vmon(slot, channel)
            now = time.time(); dt = now - last; last = now
            vset, state = ctrl.step(time.strftime('%Y-%m-%d %H:%M:%S'), imon, dt)
            if abs(vset - caen.get_ch_v0(slot, channel)) > 0.5:
                caen.set_ch_v0(slot, channel, vset)
            print(f"  {time.strftime('%H:%M:%S')} {state:<7} "
                  f"vmon={vmon:6.1f} imon={imon:6.2f} -> vset={vset:6.1f}")
            time.sleep(poll)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description='HV soft-trip training controller + backtest on real data.')
    ap.add_argument('--csv', help='hv_monitor.csv to backtest.')
    ap.add_argument('--channel', default='1:0', help="slot:ch (default 1:0).")
    ap.add_argument('--out', default=None, help='output PNG (backtest).')
    ap.add_argument('--show', action='store_true', help='open a window.')
    ap.add_argument('--v-start', type=float, default=None,
                    help='initial commanded V (default = target).')
    # controller params
    ap.add_argument('--i-comp', type=float, default=10.0)
    ap.add_argument('--i-comp-frac', type=float, default=0.95)
    ap.add_argument('--i-safe', type=float, default=2.0)
    ap.add_argument('--v-target', type=float, default=420.0)
    ap.add_argument('--v-step-up', type=float, default=10.0)
    ap.add_argument('--v-step-down', type=float, default=25.0)
    ap.add_argument('--dwell', type=float, default=30.0)
    ap.add_argument('--recover-dwell', type=float, default=30.0)
    ap.add_argument('--backoff-after', type=float, default=6.0)
    ap.add_argument('--caen-trip-time', type=float, default=30.0)
    args = ap.parse_args()

    cfg = TrainerConfig(
        i_comp=args.i_comp, i_comp_frac=args.i_comp_frac, i_safe=args.i_safe,
        v_target=args.v_target, v_step_up=args.v_step_up,
        v_step_down=args.v_step_down, dwell=args.dwell,
        recover_dwell=args.recover_dwell, backoff_after=args.backoff_after,
        caen_trip_time=args.caen_trip_time)

    if not args.csv:
        ap.error('Provide --csv to backtest (live mode is run via run_live()).')
    out = args.out
    if out is None:
        base = os.path.dirname(os.path.abspath(args.csv))
        out = os.path.join(base, f'hv_train_backtest_{args.channel.replace(":","-")}.png')
    backtest_csv(args.csv, args.channel, cfg, out, args.show, args.v_start)


if __name__ == '__main__':
    main()
