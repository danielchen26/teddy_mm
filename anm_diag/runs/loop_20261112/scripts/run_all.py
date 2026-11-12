"""Run every audit that has no record yet, two at a time: pooled first, then transfer, then the single-head audits."""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from _paths import HERE, RUN, audits

order = {"pooled": 0, "transfer": 1, "target": 2}
todo = sorted((n for n in audits() if not (RUN / n / "anm_record.json").exists()), key=lambda n: (order[audits()[n][0]], n))
print(len(todo), "audits to run", flush=True)


def one(n):
    r = subprocess.run([sys.executable, str(HERE / "run_one.py"), n], capture_output=True, text=True)
    print((r.stdout.strip().splitlines() or ["?"])[-1] if r.returncode == 0 else f"FAILED {n}: {r.stderr.strip()[-400:]}", flush=True)


with ThreadPoolExecutor(2) as ex:
    list(ex.map(one, todo))
