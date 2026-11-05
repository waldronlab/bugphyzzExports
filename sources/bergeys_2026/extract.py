#!/usr/bin/env python3 -I
"""Extract traits from Bergey's genus abstracts with the Gemini API.

Reads the API key from pi's auth.json inside this process only; never prints it.
Usage: bergeys_extract.py <system.txt> <batches_dir> <out_dir> [concurrency]
"""
import json, sys, os, glob, time, re, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

system_txt, batches_dir, out_dir = sys.argv[1:4]
conc = int(sys.argv[4]) if len(sys.argv) > 4 else 4
os.makedirs(out_dir, exist_ok=True)
key = json.load(open(os.path.expanduser("~/.pi/agent/auth.json")))["google"]["key"]
system = open(system_txt).read()
URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-flash-latest:generateContent"

def call(batch_path):
    name = os.path.basename(batch_path)[:-5]
    out = os.path.join(out_dir, name + ".json")
    if os.path.exists(out):
        return name, "cached"
    genera = json.load(open(batch_path))
    body = {
        "system_instruction": {"parts": [{"text": system}]},
        "contents": [{"role": "user", "parts": [{"text": "Extract traits for these genera and reply with the JSON list only.\n\n" + json.dumps(genera, ensure_ascii=False)}]}],
        "generationConfig": {"temperature": 0, "response_mime_type": "application/json", "max_output_tokens": 65536},
    }
    for attempt in range(5):
        try:
            req = urllib.request.Request(URL, data=json.dumps(body).encode(), headers={"Content-Type": "application/json", "x-goog-api-key": key})
            resp = json.load(urllib.request.urlopen(req, timeout=600))
            text = resp["candidates"][0]["content"]["parts"][0]["text"]
            open(os.path.join(out_dir, name + ".raw"), "w").write(text)
            data = json.loads(re.sub(r"^```(json)?|```$", "", text.strip(), flags=re.M))
            got = {d.get("doi") for d in data}; want = {g["doi"] for g in genera}
            if want - got:
                missing = want - got
                return name, f"parsed {len(data)} but missing {len(missing)} dois"
            json.dump(data, open(out, "w"), ensure_ascii=False)
            return name, f"ok {len(data)}"
        except (urllib.error.HTTPError, urllib.error.URLError, json.JSONDecodeError, KeyError, IndexError) as e:
            msg = getattr(e, "code", None) or type(e).__name__
            time.sleep(10 * (attempt + 1))
            last = f"{msg}: {str(e)[:120]}"
    return name, "FAILED " + last

batches = sorted(glob.glob(os.path.join(batches_dir, "batch_*.json")))
with ThreadPoolExecutor(conc) as ex:
    futs = {ex.submit(call, b): b for b in batches}
    for f in as_completed(futs):
        name, status = f.result()
        print(f"{name}: {status}", flush=True)
