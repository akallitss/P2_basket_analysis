# det2 no-hit checks (27 Nov 2026), run from cosmic_bench_analysis/ after 34 and 39

| script | question | result on det2_long_repos |
|---|---|---|
| m3tests.py | are no-hit tracks fakes / kinked / pile-up? (top vs bottom M3 pair slopes, chi2, NClus, dt to neighbours) | no: identical to hit tracks at every angle |
| spots.py, spots2.py | do no-hit tracks cluster in fixed spots? | only at 4 big pillars (4 % of no-hits); vertical no-hits are spatially random |
| masktest.py | no-hit vs distance to pillars of Saclay V1/V2 and CERN masks | no shadow for any design |
| offscan.py, offscan2.py | shift the Gerber pillar pattern by (dx,dy) over one lattice cell, globally and per quadrant | best z 2.6-4.7 of ~1600 trials, inconsistent offsets: no shadow |
| strokes.py | what are the non-pillar strokes in the masks? | outlines, frame, text only |
