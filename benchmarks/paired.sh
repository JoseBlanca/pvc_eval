# The wall time and peak RSS of str-regions from two worktrees, run alternately
# (benchmarks/README.md, "The STR finder").
# Usage: bash benchmarks/paired.sh WORKTREE_A WORKTREE_B INPUT_FASTA PAIRS [STR_REGIONS_OPTIONS...]
# Runs str-regions from A and B alternately, PAIRS times, A first in odd pairs and B first in
# even ones, so that a drift of the machine favours each side as often. It prints each side's
# wall time, user+sys time, peak RSS and the 1-minute load average before the run, then the
# median ratio B/A of the wall and of the CPU times. The two outputs of every pair are compared.
A=$1; B=$2; IN=$3; N=$4; shift 4
O=.profiles/paired
mkdir -p $O
run() {
  local project=$1 side=$2; shift 2
  load=$(sysctl -n vm.loadavg | awk '{print $2}')
  /usr/bin/time -l uv run --no-sync --project $project pvc_eval str-regions --reference $IN --output $O/$side.bed "$@" 2>&1 \
    | awk -v load=$load '/real/{t=$1; c=$3+$5} /maximum resident/{m=$1} END{print t, c, m, load}'
}
: > $O/times
for i in $(seq $N); do
  rm -f $O/a.bed $O/b.bed
  if (( i % 2 )); then a=$(run $A a "$@"); b=$(run $B b "$@"); else b=$(run $B b "$@"); a=$(run $A a "$@"); fi
  [[ -s $O/a.bed ]] || { echo "no output from A"; exit 1; }
  cmp -s $O/a.bed $O/b.bed && same=equal || same=DIFFERENT
  echo "$i A $a B $b $same" | tee -a $O/times
done
python3 - $O/times <<'PY'
import sys, statistics
rows = [l.split() for l in open(sys.argv[1])]
a, ca, ra = ([float(r[k]) for r in rows] for k in (2, 3, 4))
b, cb, rb = ([float(r[k]) for r in rows] for k in (7, 8, 9))
loads = [float(r[k]) for r in rows for k in (5, 10)]
sd = lambda xs: statistics.stdev(xs) if len(xs) > 1 else 0.0
print(f"A wall {statistics.mean(a):.3f} s sd {sd(a):.3f}, cpu {statistics.mean(ca):.3f} s"
      f"  B wall {statistics.mean(b):.3f} s sd {sd(b):.3f}, cpu {statistics.mean(cb):.3f} s")
ratios = [y / x for x, y in zip(a, b)]
print(f"median B/A wall {statistics.median(ratios):.4f}  min {min(ratios):.4f} max {max(ratios):.4f}"
      f"  cpu {statistics.median(y / x for x, y in zip(ca, cb)):.4f}")
print(f"peak RSS median A {statistics.median(ra)/1e6:.1f} MB  B {statistics.median(rb)/1e6:.1f} MB"
      f"  load {min(loads):.1f} to {max(loads):.1f}"
      f"  outputs {'all equal' if all(r[-1] == 'equal' for r in rows) else 'DIFFERENT'}")
PY
