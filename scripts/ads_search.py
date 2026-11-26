"""NASA ADS full-text search for a list of TOIs (background check of flagged candidates).

    python scripts/ads_search.py OUT.csv TOI [TOI ...]

The API token is read from ~/.ads/dev_key (outside the repo; never printed or committed).
For each TOI (e.g. 4842.01) it searches full text + abstract + title for:
  "TOI-4842", "TOI 4842", "TOI4842", "TIC 131034657", "TIC131034657"
(the TIC ID comes from results/hunt_v2_2026-11-25.csv). Returns bibcode, title, year, and the
highlighted text snippets ADS returns, so you can see what each paper says.
"""
import csv, json, os, sys, time, urllib.parse, urllib.request
import pandas as pd

KEY = open(os.path.expanduser("~/.ads/dev_key")).read().strip()
API = "https://api.adsabs.harvard.edu/v1/search/query"


def query(q):
    params = {"q": q, "fl": "id,bibcode,title,year,doctype,pub", "rows": 50,
              "hl": "true", "hl.fl": "body,abstract,title", "hl.snippets": 3, "hl.fragsize": 200}
    req = urllib.request.Request(API + "?" + urllib.parse.urlencode(params),
                                 headers={"Authorization": f"Bearer {KEY}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


if __name__ == "__main__":
    out, tois = sys.argv[1], [float(x) for x in sys.argv[2:]]
    hunt = pd.read_csv("results/hunt_v2_2026-11-25.csv")
    rows = []
    for toi in tois:
        n = str(int(toi))
        tic = str(hunt[hunt.TOI.round(2) == round(toi, 2)].star.iloc[0]).split()[-1]
        terms = [f'"TOI-{n}"', f'"TOI {n}"', f'"TOI{n}"', f'"TIC {tic}"', f'"TIC{tic}"']
        q = " OR ".join(f"full:{t}" for t in terms)
        d = query(q)
        docs, hl = d["response"]["docs"], d.get("highlighting", {})
        print(f"TOI {toi:.2f} (TIC {tic}): {d['response']['numFound']} papers", flush=True)
        for doc in docs:
            h = hl.get(doc.get("id", ""), {}) or {}
            snip = " … ".join(sum(h.values(), []))[:600] if h else ""
            rows.append(dict(TOI=toi, tic=tic, bibcode=doc["bibcode"], year=doc.get("year"),
                             title=(doc.get("title") or [""])[0], pub=doc.get("pub"),
                             doctype=doc.get("doctype"), snippet=snip))
            print(f"   {doc['bibcode']}  {(doc.get('title') or [''])[0][:90]}", flush=True)
        if not docs:
            rows.append(dict(TOI=toi, tic=tic, bibcode="", title="(no papers found)"))
        time.sleep(1)
    pd.DataFrame(rows).to_csv(out, index=False)
