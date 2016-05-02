#!/home/maxhutch/anaconda3/bin/python3

import json
from os import getcwd
from os.path import join
from nekpy.dask.subgraph import series
from nekpy.dask.tasks import configure
from nekpy.dask.utils import outer_product, work_name
from nekpy.dask.tasks import metal
from nekpy.dask import run_all
from dask.multiprocessing import get
from dask.async import get_sync
from dask import set_options
from importlib import import_module

metal = import_module(".metal_slurm", "nekpy.dask")


from sys import argv
with open(argv[1], "r") as f:
    base = json.load(f)

with open(argv[2], "r") as f:
    sweeps = json.load(f)

with open(argv[3], "r") as f:
    tusr = f.read()

base["prefix"] = sweeps["prefix"]
del sweeps["prefix"]

# Take simple outer product of contents of sweep file
candidates = list(outer_product(sweeps))

# Filter out the cases we don't want
overrides = []
for c in candidates:
    if c["order"] * c["elms"] > 256:
        continue
    if c["order"] * c["elms"] < 32:
        continue
    overrides.append(c)

# Tune the remaining cases
for ov in overrides:
    ov["name"] = work_name(base["prefix"], ov)
    ov["shape_mesh"] = [ov["elms"], ov["elms"], 4*ov["elms"]]
    nodes = max(1, int(4 * (ov["elms"])**3 / 32))
    ov["procs"] = 32*nodes
    ov["io_files"] = -nodes
    ov["dt"] = (2/(ov["elms"]*(ov["order"]-1)**2))/0.0558519

workdirs = [join(getcwd(), x["name"]) for x in overrides]
configs = [configure(base, override, workdir) for override, workdir in zip(overrides, workdirs)]
res = [series(config, tusr, job_time = 8.0) for config in configs]
run_all(res, base, get=get, num_workers=2)

