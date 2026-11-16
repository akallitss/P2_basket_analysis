#!/usr/bin/env python3
"""Closure tests of channel_efficiency.py (run: python3 test_channel_efficiency.py, or pytest). Simulated muons with a
known efficiency must come back within the quoted error whatever the tracking resolution, also when the predicted
points sit on the pad lattice; random hits must be subtracted; resolution tails must show up in the R_lo scan; impure
tags must NOT be corrected (they bias by their fraction: why only scintillator-confirmed tags are used)."""
import channel_efficiency as CE


def _check(r, eps, nsig=3.0):
    assert abs(r['eps_est'] - eps) < nsig * r['err'], r


def test_resolution_is_removed():
    for sg in (3.0, 8.0, 12.0):
        r = CE.toy(sg, 0.85, 12000, seed=3)
        _check(r, 0.85)
        assert r['naive'] < 0.6          # the naive own-pad number is far off, as in the data


def test_lattice_predictions():
    for sg in (3.0, 8.0):
        _check(CE.toy(sg, 0.85, 12000, seed=4, quantize=True), 0.85)


def test_random_hits_are_subtracted():
    _check(CE.toy(8.0, 0.7, 12000, noise=0.0005, seed=5), 0.7)


def test_tails_show_in_the_rlo_scan():
    r = CE.toy(3.0, 0.85, 12000, tail=(0.15, 20.0), seed=6)
    assert r['eps_rlo56'] - r['eps_rlo24'] > 0.03, r     # the scan rises: not all muons inside at small R_lo
    assert abs(r['eps_rlo56'] - 0.85) < 3 * r['err'], r  # and recovers at large R_lo


def test_impure_tags_bias_by_the_impurity():
    r = CE.toy(3.0, 0.85, 12000, impure=0.3, seed=7)
    assert abs(r['eps_est'] - 0.85 * 0.7) < 3 * r['err'], r


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print('ok', name)
