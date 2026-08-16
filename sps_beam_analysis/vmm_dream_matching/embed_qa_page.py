#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
embed_qa_page.py -- put the QA data inside the static QA page.

The website is static HTML with no server and no CDN, so the page carries its
own data: this replaces the contents of

    <script type="application/json" id="qa-data"> ... </script>

with the JSON that build_qa_data.py produced. Idempotent -- the page is both
the template and the output, so run it again after every rebuild.

    python3 build_qa_data.py --out sps_qa_data.json          # on lxplus
    python3 embed_qa_page.py sps_qa_data.json \\
            ~/Documents/PostDocSaclay/cern-site/notes/2026-08-16-sps-qa.html
"""
import argparse
import datetime as dt
import json
import os
import re
import sys

OPEN = '<script type="application/json" id="qa-data">'
CLOSE = "</script>"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("page")
    ap.add_argument("--stamp", action="store_true", default=True,
                    help="refresh the <meta name=updated> date")
    args = ap.parse_args()

    with open(args.data) as f:
        payload = f.read().strip()
    json.loads(payload)                      # fail loudly on a broken build
    # </script> inside the data would end the block early; JSON escapes it
    payload = payload.replace("</", "<\\/")

    with open(args.page) as f:
        html = f.read()
    i = html.find(OPEN)
    if i < 0:
        sys.exit(f"{args.page}: no qa-data script block to fill")
    j = html.find(CLOSE, i)
    html = html[:i + len(OPEN)] + payload + html[j:]

    if args.stamp:
        today = dt.date.today().isoformat()
        html = re.sub(r'(<meta name="updated" content=")[^"]*(">)',
                      lambda m: m.group(1) + today + m.group(2), html)

    with open(args.page, "w") as f:
        f.write(html)
    print(f"embedded {len(payload)/1e6:.2f} MB into {args.page} "
          f"({os.path.getsize(args.page)/1e6:.2f} MB total)")


if __name__ == "__main__":
    main()
