# Parallelism

**Purpose.** pvc_eval runs on the free-threaded Python 3.15t, which the
owner chose on 19 October 2026 so that its threads run Python code at the
same time. This category says where threads would pay, what they cost, and,
once there are threads, whether they work as well as they could.

**Sent.** Always once pvc_eval uses threads. Before that, when the profile
shows independent work that takes most of the time.

## Rules

- **Threads, not processes.** On 3.15t, threads share the reference, the
  records and the code, with no copy and no pickling, and run on all the
  cores. Processes are kept for what threads cannot do, and are the
  owner's to decide: each is a whole interpreter with its own copy of what
  it needs.
- **Where it splits.** The stages are a stream in order, so the natural
  unit is the contig: contigs share no locus. A finding says how the work
  is spread over contigs on the workload, since one large contig among
  small ones leaves the other threads idle at the end, and how the output
  is put back in the order of the `.fai`: the table must not depend on
  which thread finished first. `ThreadPoolExecutor.map` gives the results
  in the order of the input, which is the order to want.
- **Shared state.** What threads share is read only: the frozen records,
  the reference, the settings. Anything a thread writes is its own until
  it returns it. A finding that proposes threads names every object two
  threads would touch and says why each is safe, and goes to the code
  review for that call.
- **A lock in the hot path is a finding.** A lock taken per record or per
  locus makes the threads take turns. Tachyon over all threads in `wall`
  and `cpu` mode shows it: threads whose wall time is much larger than
  their CPU time are waiting on something.
- **Scaling.** The measurement of a change that adds threads, or of a
  stage that has them, is the wall time and the peak memory at one, two,
  four, eight and sixteen threads, the output compared with the run on one
  thread byte for byte. Speed that stops growing well before the number
  of cores means a part that runs alone, a lock, or work spread unevenly,
  and the profile says which.
- **Memory grows with the threads.** Each thread holds what its task
  holds, so the peak at eight threads is about eight times the peak of one
  task, plus what is shared. Both peaks are stated.
- **Something that turns the GIL back on.** A C extension imported on
  3.15t that has not declared that it runs without the GIL turns the GIL
  on for the whole process, and the threads run one at a time with no
  error, only a `RuntimeWarning`. A run that does not scale is checked for
  it first: `sys._is_gil_enabled()` after the imports.
- **What it costs in code.** A pool, the order of the output, a test that
  runs with one thread and with several and compares the outputs as equal,
  errors raised inside a thread brought back to the caller. Said in full
  under the effect on the code, and judged by the table of `SKILL.md`.
