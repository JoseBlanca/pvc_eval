"""Where the peak memory of `pvc_eval str-regions` comes from: the peak bytes tracemalloc
traced, and the source lines that held the most when the bytes held were near their
highest.

    uv run --no-sync python benchmarks/peak_str_regions.py --reference FASTA --output BED

The arguments are those of `pvc_eval str-regions`. A second thread takes a snapshot
whenever the bytes held have grown by _SNAPSHOT_GROWTH over the last one, so the snapshot
printed was taken within that share of the highest point it saw. tracemalloc slows the
run, so its time is not a measurement.
"""

import sys
import threading
import tracemalloc

from pvc_eval.cli import main

_NUM_TOP_LINES = 15
_NUM_FRAMES = 1
_SNAPSHOT_GROWTH = 1.05
_POLL_SECONDS = 0.05
_MIB = 1024 * 1024


def _watch(done: threading.Event, largest: list[tuple[int, tracemalloc.Snapshot]]) -> None:
    """Keep in largest the snapshot taken when the bytes held were highest."""
    while not done.wait(_POLL_SECONDS):
        held, _ = tracemalloc.get_traced_memory()
        if not largest or held > largest[0][0] * _SNAPSHOT_GROWTH:
            largest[:] = [(held, tracemalloc.take_snapshot())]


def run() -> None:
    tracemalloc.start(_NUM_FRAMES)
    done = threading.Event()
    largest: list[tuple[int, tracemalloc.Snapshot]] = []
    watcher = threading.Thread(target=_watch, args=(done, largest))
    watcher.start()
    try:
        main(["str-regions", *sys.argv[1:]])
    finally:
        done.set()
        watcher.join()
    _, peak = tracemalloc.get_traced_memory()
    print(f"peak traced: {peak / _MIB:.1f} MiB")
    held, snapshot = largest[0]
    print(f"held when the snapshot was taken: {held / _MIB:.1f} MiB")
    print(f"top {_NUM_TOP_LINES} lines then:")
    for statistic in snapshot.statistics("lineno")[:_NUM_TOP_LINES]:
        print(statistic)


if __name__ == "__main__":
    run()
