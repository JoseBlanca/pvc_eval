# Where pvc_eval's STR finder comes from

*19 October 2026.* pvc_eval's STR finder (`doc/spec.md` §5.1, `doc/specs/str_finder.md`) is a
port of pvc's, from the repository `pop_var_caller` at commit
`bbcc9f4658d57a2bbedc97b0edfcaff366205bb3`, of 18 October 2026, its HEAD on the day the port began.
pvc_eval imports nothing from pvc. A change to the port that follows a later commit of pvc updates
this record in the same commit.

| step of the port | pvc file | pvc function or constant |
|---|---|---|
| 1. candidates | `src/tandem_repeat.rs` | `find_tandem_repeats`, `maximal_scoring_subsequences`, `is_primitive_motif`, `CANONICAL`, `bases_match`; `DEFAULT_MATCH_REWARD`, `DEFAULT_MISMATCH_PENALTY` |
| 1. copies a candidate needs | `src/cli/repeat_catalog.rs` | `SCANNER_MIN_COPIES` |
| 2. contig-end flank | `src/repeat_catalog/row.rs`, `src/repeat_catalog/criteria.rs` | `row_for_interval`; `CATALOG_MIN_FLANK_BP` |
| 3. copy floor on the detected span | `src/region_typing/segment_criteria.rs` | `prefilter`; pvc_eval's default floors are `CATALOG_MIN_COPIES` of `src/repeat_catalog/criteria.rs` |
| 4. scope tests | `src/region_typing/segment_criteria.rs` | `classify`, its loop before `split_bundles`; `is_compound`, `upper` |
| 5. bundles | `src/region_typing/segment_criteria.rs` | `split_bundles`, `joins_cluster`; `DEFAULT_BUNDLE_THRESHOLD` |
| 6. STR segments | `src/region_typing/segment_criteria.rs` | `finish_locus`, `minimal_trim`, `recompute_purity`, `bundle_member_from_bases`; `DEFAULT_MIN_PURITY`, `DEFAULT_MIN_SCORE`, `DEFAULT_MIN_PERIOD`, `DEFAULT_MAX_PERIOD` |
| 7. satellites and bundles | `src/region_typing/mod.rs` | `resolve_features`, `absorb_into` (the fixed point that `absorb_by_iteration` runs in pvc's tests, not `absorb_by_sweep`), `merge_runs`, `coverage_runs`; `bundle_clusters` of `segment_criteria.rs`; `DEFAULT_MAX_STR_LEN`, `DEFAULT_MAX_STR_BUNDLE_LEN` |
| 8. generic regions | `src/region_typing/mod.rs` | `fill_generic_gaps` |
| the whole contig | `src/region_typing/mod.rs` | `partition_resident_in`, `partition_resident` |
| `str-regions --regions` | `src/region_typing/mod.rs` | `emit_into`, `clips_at_a_bed_edge` |
| tests and fixtures | `src/tandem_repeat.rs`, `src/region_typing/segment_criteria.rs`, `src/region_typing/mod.rs`, `src/repeat_catalog/row.rs`; `tests/data/tandem_repeat/` | the tests listed in `doc/specs/str_finder_tests.md`; `synthetic_ref.fa`, `golden.ssr_catalog.bed.gz` |
