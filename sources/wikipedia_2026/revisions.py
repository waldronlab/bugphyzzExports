"""Find the English Wikipedia revision each introduction was taken from.

taxoniq fetched the introductions live, between WINDOW_START and WINDOW_END.
For each title, the revision current at WINDOW_END is looked up; if the
article was not edited inside the window, that revision is exactly the text
that was used. For articles edited inside the window, the introduction of
each revision in the window is compared with the text taxoniq stored, and
the closest one is kept (match_ratio; matched_by = "text").

Usage: python3 -I revisions.py <records.jsonl> <out.tsv>
"""

import csv
import difflib
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request

WINDOW_START = "2026-10-26T00:00:00Z"
WINDOW_END = "2026-10-30T16:05:00Z"
API = "https://en.wikipedia.org/w/api.php"
UA = "bugphyzz-wikipedia-revisions/1.0 (https://github.com/waldronlab/bugphyzzExports)"


def get(params):
    params = dict(params, format="json", formatversion="2")
    url = API + "?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as res:
                return json.load(res)
        except Exception as e:  # retry network errors
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))


def latest(titles):
    """Latest revision for up to 50 titles; follows redirects."""
    doc = get({"action": "query", "prop": "revisions", "rvprop": "ids|timestamp",
               "redirects": "1", "titles": "|".join(titles)})
    q = doc["query"]
    alias = {}
    for k in ("normalized", "redirects"):
        for x in q.get(k, []):
            alias[x["from"]] = x["to"]
    pages = {p["title"]: p for p in q["pages"]}
    out = {}
    for t in titles:
        r = t
        while r in alias:
            r = alias[r]
        out[t] = (r, pages.get(r))
    return out


def at(title, when):
    """Revision of one page current at a timestamp."""
    doc = get({"action": "query", "prop": "revisions", "rvprop": "ids|timestamp",
               "rvlimit": "1", "rvdir": "older", "rvstart": when, "titles": title})
    revs = doc["query"]["pages"][0].get("revisions", [])
    return revs[0] if revs else None


def in_window(title):
    """Revisions current at any time in the window, oldest first."""
    doc = get({"action": "query", "prop": "revisions", "rvprop": "ids|timestamp",
               "rvlimit": "50", "rvdir": "newer", "rvstart": WINDOW_START,
               "rvend": WINDOW_END, "titles": title})
    revs = doc["query"]["pages"][0].get("revisions", [])
    first = at(title, WINDOW_START)
    return ([first] if first else []) + revs


def norm(text):
    for tag in ("style", "sup", "table"):
        text = re.sub(rf"<{tag}\b.*?</{tag}>", " ", text, flags=re.S)
    text = html.unescape(re.sub(r"<[^>]+>", "", text))
    text = re.sub(r"\[\d+\]|\[citation needed\]", " ", text)
    text = re.sub(r"\s+([,.;:)])", r"\1", text)
    return re.sub(r"\s+", " ", text).strip()


def intro(revid):
    doc = get({"action": "parse", "oldid": revid, "prop": "text", "section": "0"})
    raw = doc["parse"]["text"]
    paras = re.findall(r"<p[^>]*>(.*?)</p>", raw, flags=re.S)
    return norm(" ".join(paras))


def closest(title, text):
    """(revision, ratio, how): the revision whose introduction is closest to
    the stored text. If the introduction did not change in the window, the
    revision current at the end is kept."""
    target = norm(text)[:1500]
    scored = []
    for rev in in_window(title):
        it = intro(rev["revid"])[:len(target) + 200]
        scored.append((rev, difflib.SequenceMatcher(None, target, it).ratio(), it))
        time.sleep(0.1)
    if not scored:
        return None
    if len({it for _, _, it in scored}) == 1:
        rev, ratio, _ = scored[-1]
        return rev, ratio, "intro unchanged"
    top = max(r for _, r, _ in scored)
    rev, ratio, _ = [x for x in scored if x[1] == top][-1]
    return rev, ratio, "text"


def main(records_path, out_path):
    recs = [json.loads(line) for line in open(records_path)]
    titles = [r["en_wiki_title"] for r in recs]
    found = {}
    for i in range(0, len(titles), 50):
        found.update(latest(titles[i:i + 50]))
        time.sleep(0.2)
    rows = []
    n_single = 0
    for r in recs:
        resolved, page = found[r["en_wiki_title"]]
        row = {"taxid": r["taxid"], "en_wiki_title": r["en_wiki_title"],
               "resolved_title": resolved, "pageid": "", "revid": "",
               "rev_timestamp": "", "revid_at_window_start": "",
               "edited_in_window": "", "matched_by": "", "match_ratio": "",
               "note": ""}
        if page is None or page.get("missing") or not page.get("revisions"):
            row["note"] = "page not found"
            rows.append(row)
            continue
        row["pageid"] = page["pageid"]
        rev = page["revisions"][0]
        if rev["timestamp"] > WINDOW_START:
            # edited since the fetch began: look the revisions up at both ends
            n_single += 1
            rev = at(resolved, WINDOW_END)
            start = at(resolved, WINDOW_START)
            time.sleep(0.1)
            if rev is None:
                row["note"] = "page created after the fetch"
                rows.append(row)
                continue
            row["revid_at_window_start"] = start["revid"] if start else ""
            edited = start is None or start["revid"] != rev["revid"]
            row["edited_in_window"] = "TRUE" if edited else "FALSE"
            row["matched_by"] = "date"
            if edited:
                best = closest(resolved, r["text"])
                if best:
                    rev = best[0]
                    row["matched_by"] = best[2]
                    row["match_ratio"] = f"{best[1]:.3f}"
        else:
            row["edited_in_window"] = "FALSE"
            row["matched_by"] = "date"
        row["revid"] = rev["revid"]
        row["rev_timestamp"] = rev["timestamp"]
        rows.append(row)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n_edit = sum(r["edited_in_window"] == "TRUE" for r in rows)
    n_miss = sum(bool(r["note"]) for r in rows)
    by_text = [float(r["match_ratio"]) for r in rows if r["matched_by"] == "text"]
    print(f"{len(rows)} records; {n_single} looked up by date; "
          f"{n_edit} edited in the window, of which {len(by_text)} changed the "
          f"introduction (lowest match_ratio {min(by_text, default=1):.3f}); {n_miss} without a revision")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
