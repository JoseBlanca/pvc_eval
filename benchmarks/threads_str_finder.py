"""How segment_contig scales across threads (benchmarks/README.md, "The STR finder").

    uv run --no-sync python benchmarks/threads_str_finder.py FASTA THREADS...

For each number of threads N, N threads each segment the first contig of FASTA at once,
after a barrier, and the script prints the wall time of each call and the CPU time of the
process divided by N times that of one call on one thread. Contention between the threads,
such as the shared closure cell of a frozen dataclass's __init__ on 3.15t, shows as that
ratio growing with N; the load of other processes moves the wall times but hardly the ratio.
"""

import argparse
import sys
import threading
import time
from pathlib import Path

from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import segment_contig


def _time_threads(
    reference: IndexedFasta, contig: str, num_threads: int
) -> tuple[list[float], float]:
    barrier = threading.Barrier(num_threads)
    walls = [0.0] * num_threads

    def segment(index: int) -> None:
        barrier.wait()
        started = time.perf_counter()
        for _ in segment_contig(reference, contig):
            pass
        walls[index] = time.perf_counter() - started

    threads = [threading.Thread(target=segment, args=(index,)) for index in range(num_threads)]
    cpu_started = time.process_time()
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    return walls, time.process_time() - cpu_started


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("threads", nargs="+", type=int)
    arguments = parser.parse_args()
    with IndexedFasta(arguments.fasta) as reference:
        contig = reference.contigs[0]
        _time_threads(reference, contig, 1)  # warm up
        _, one_cpu = _time_threads(reference, contig, 1)
        for num_threads in arguments.threads:
            walls, cpu = _time_threads(reference, contig, num_threads)
            print(
                f"{num_threads} threads: {min(walls):.2f} to {max(walls):.2f} s per call,"
                f" cpu {cpu:.2f} s, {cpu / (num_threads * one_cpu):.2f} x one thread's per call"
                f" (one thread {one_cpu:.2f} s cpu), GIL {'on' if sys._is_gil_enabled() else 'off'}",
                flush=True,
            )


if __name__ == "__main__":
    main()
