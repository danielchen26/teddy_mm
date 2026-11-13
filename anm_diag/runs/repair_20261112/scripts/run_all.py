"""Run every audit of one split that has no record yet, two at a time, the pooled audits first:  python run_all.py <seed>"""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from _paths import HERE, RUN, audits

seed = int(sys.argv[1])
todo = sorted((n for n, a in audits().items() if a[3] == seed and not (RUN / n / "anm_record.json").exists()), key=lambda n: (audits()[n][0] != "pooled", n))
print(len(todo), "audits to run", flush=True)


def one(n):
    r = subprocess.run([sys.executable, str(HERE / "run_one.py"), n], capture_output=True, text=True)
    print((r.stdout.strip().splitlines() or ["?"])[-1] if r.returncode == 0 else f"FAILED {n}: {r.stderr.strip()[-400:]}", flush=True)


with ThreadPoolExecutor(2) as ex:
    list(ex.map(one, todo))
