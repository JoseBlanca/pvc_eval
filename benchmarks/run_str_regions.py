"""Run `pvc_eval str-regions` with the arguments given, as the console script does, so
that a profiler or tracemalloc can run it as a script (pvc_eval has no __main__).

    uv run --no-sync python benchmarks/run_str_regions.py --reference FASTA --output BED
"""

from pvc_eval.cli import main

if __name__ == "__main__":
    main(["str-regions", *__import__("sys").argv[1:]])
