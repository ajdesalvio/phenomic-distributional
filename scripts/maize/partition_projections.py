"""Stream the original projection CSV into disposable per-heldout caches.

Byte-preserving rows, SHA256 source provenance, no model or input modification.
Rerunning with an existing complete manifest validates content hashes and size/mtime before reuse.
"""
import csv
import hashlib
import json
import pathlib
import sys
import time

def file_hash(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

source = pathlib.Path(sys.argv[1]).resolve()
out = pathlib.Path(sys.argv[2]).resolve()
manifest = out / "partition_manifest.json"
stat = source.stat()
if manifest.exists():
    prior = json.loads(manifest.read_text())
    assert prior["source"] == str(source)
    assert prior["source_size"] == stat.st_size
    assert prior["source_sha256"] == file_hash(source)
    assert all(file_hash(out / x["file"]) == x["sha256"] for x in prior["partitions"])
    assert prior["source_mtime_ns"] == stat.st_mtime_ns
    assert all((out / x["file"]).stat().st_size == x["bytes"] for x in prior["partitions"])
    print("Reusing complete content-verified partition cache.", flush=True)
    sys.exit(0)
if out.exists() and any(out.iterdir()):
    raise RuntimeError("Incomplete/nonempty cache: choose a new output directory.")
out.mkdir(parents=True, exist_ok=True)
handles, counts = {}, {}
sha = hashlib.sha256()
started = time.time()
try:
    with source.open("rb", buffering=4 * 1024 * 1024) as src:
        header = src.readline()
        sha.update(header)
        assert header.split(b",")[6] == b"Heldout.Env"
        for n, line in enumerate(src, 1):
            sha.update(line)
            first = line.split(b",", 7)
            env = (next(csv.reader([line.decode("utf-8")]))[6]
                   if b'"' in b",".join(first[:7]) else first[6].decode("utf-8"))
            if not env or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-" for c in env):
                raise ValueError(f"Unsafe heldout identifier: {env!r}")
            if env not in handles:
                handles[env] = (out / f"{env}.csv").open("wb", buffering=4 * 1024 * 1024)
                handles[env].write(header)
                counts[env] = 0
            handles[env].write(line)
            counts[env] += 1
            if n % 1000000 == 0:
                print(f"Partitioned {n:,} rows in {time.time()-started:.0f}s", flush=True)
finally:
    for h in handles.values():
        h.close()
result = dict(source=str(source), source_size=stat.st_size,
              source_mtime_ns=stat.st_mtime_ns, source_sha256=sha.hexdigest(),
              partitions=[dict(heldout=e, file=f"{e}.csv", rows=counts[e],
                               bytes=(out / f"{e}.csv").stat().st_size, sha256=file_hash(out / f"{e}.csv"))
                          for e in sorted(counts)])
manifest.write_text(json.dumps(result, indent=2))
print(f"Completed {sum(counts.values()):,} rows across {len(counts)} environments.", flush=True)
