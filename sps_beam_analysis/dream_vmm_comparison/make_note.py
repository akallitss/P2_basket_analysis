#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_note.py -- wrap a report.html fragment into a self-contained note.

The site's notes must carry no external references (they are precached for
offline reading), so every <img src="figures/x.png"> is inlined as a data: URI.

    python3 make_note.py report.html report_note.html \
        --title "..." --summary "..." --tags "a, b, c"
"""
import argparse
import base64
import os
import re

ap = argparse.ArgumentParser()
ap.add_argument("src")
ap.add_argument("dst")
ap.add_argument("--title", required=True)
ap.add_argument("--summary", required=True)
ap.add_argument("--tags", required=True)
a = ap.parse_args()

HERE = os.path.dirname(os.path.abspath(a.src)) or "."
body = open(a.src).read()

# the fragment already carries <title>, <meta name=description> and <style>
title = re.search(r"<title>(.*?)</title>", body, re.S)
desc = re.search(r'<meta name="description" content="(.*?)">', body, re.S)


def inline(m):
    p = os.path.join(HERE, m.group(1))
    with open(p, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode()
    return f'src="data:image/png;base64,{b64}"'


body = re.sub(r'src="([^"]+\.png)"', inline, body)

head = f"""<!--note
title: {a.title}
tags: {a.tags}
summary: {a.summary}
-->
<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
"""
out = head + body + "\n</body></html>\n"
# the fragment's <title>/<meta>/<style> belong in <head>; close it before <main>
out = out.replace("<main>", "</head><body>\n<main>", 1)
open(a.dst, "w").write(out)
print(f"wrote {a.dst}  ({os.path.getsize(a.dst) / 1024:.0f} kB)"
      f"  title={title.group(1) if title else '?'}")
