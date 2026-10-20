# The wall time and peak RSS of str-regions from two worktrees, run alternately
# (benchmarks/README.md, "The STR finder").
# Usage: bash benchmarks/paired.sh WORKTREE_A WORKTREE_B INPUT_FASTA PAIRS
# Runs str-regions from A and B alternately, PAIRS times, and prints each wall time, peak RSS,
# and the median ratio B/A. The outputs of the last pair are compared with cmp.
A=$1; B=$2; IN=$3; N=$4
O=.profiles/paired
mkdir -p $O
run() { /usr/bin/time -l uv run --no-sync --project $1 pvc_eval str-regions --reference $IN --output $O/$2.bed 2>&1 | awk '/real/{t=$1} /maximum resident/{m=$1} END{print t, m}'; }
: > $O/times
for i in $(seq $N); do
  a=$(run $A a); b=$(run $B b)
  echo "$i A $a B $b" | tee -a $O/times
done
cmp $O/a.bed $O/b.bed && echo "outputs equal"
python3 - $O/times <<'PY'
import sys, statistics
rows=[l.split() for l in open(sys.argv[1])]
a=[float(r[2]) for r in rows]; b=[float(r[5]) for r in rows]
ra=[float(r[3]) for r in rows]; rb=[float(r[6]) for r in rows]
print(f"A mean {statistics.mean(a):.3f} s sd {statistics.stdev(a):.3f}  B mean {statistics.mean(b):.3f} s sd {statistics.stdev(b):.3f}")
ratios=[y/x for x,y in zip(a,b)]
print(f"median B/A {statistics.median(ratios):.4f}  min {min(ratios):.4f} max {max(ratios):.4f}")
print(f"peak RSS median A {statistics.median(ra)/1e6:.1f} MB  B {statistics.median(rb)/1e6:.1f} MB")
PY
