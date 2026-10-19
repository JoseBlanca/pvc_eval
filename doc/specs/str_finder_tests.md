# pvc's tests of the STR finder, and what the port does with each

*19 October 2026. The list behind "How it is verified" in `doc/specs/str_finder.md`, for whoever
writes the tests. pvc is at bbcc9f4 (`doc/port_provenance.md`); line numbers are of that commit.*

Each test of pvc that the port carries is rewritten in Python with the same inputs and the same
expected outputs. Its inline sequences are copied as they are. A 1-based inclusive position of pvc,
`[s, e]`, becomes `[s - 1, e)`. A test that drops an assertion says which one in a comment, with
the reason given here. The step functions are private functions of `pvc_eval.segments`, and the
tests call them by the names their steps have in the spec.

Three configurations recur:

- **params**: pvc's `params()` (segment_criteria.rs 1516), `matched_params(0.8, 0, 5)`. Periods 2
  to 6, copy floors (10, 5, 4, 3, 3, 3), purity 0.8, score floor 0, bundle distance 5.
- **catalog**: pvc's `catalog_config()` (mod.rs 1824). Periods 2 to 6, copy floors
  (10, 5, 4, 3, 3, 3), bundle distance 50, purity 0.8, score floor 0, satellite cap 1000, bundle
  cap unbounded, scan 2/7/2.
- **pvc-default**: what pvc's `Default` gives, which pvc_eval's defaults do not. Periods 1 to 6,
  copy floors (8, 6, 6, 6, 5, 4), purity 0.8, score floor 0, bundle distance 15, satellite cap 100,
  bundle cap 100, scan 2/7/2. A test marked with it passes these values explicitly.

Every test of pvc's partition runs with a contig-end flank of 0, which is pvc's scanning route.
The contig-end flank is a step of its own (step 2), not part of `classify`'s port. So the tests of
`classify` with 5-base flanks are unchanged by it.

Class **A** is ported as is. Class **B** is ported with a change, given in the note. Class **C**
is not ported (**Open 1** of the spec).

## tandem_repeat.rs, 35 tests: 17 ported

| test | line | class | note |
|---|---|---|---|
| period_range_accepts_valid_ranges | 1147 | B | (1,6), (6,6) and (2,6) make valid `StrFinderParameters` |
| period_range_rejects_zero_min | 1164 | B | (0,6) and (0,0) raise ValueError |
| period_range_rejects_min_above_max | 1171 | B | (3,2) raises ValueError |
| default_periods_is_one_to_six | 1181 | A | |
| scan_params_default_is_two_seven_two | 1188 | B | the defaults of `StrFinderParameters` are 2, 7, 2 |
| ruzzo_tompa_matches_brute_force_on_crafted_cases | 1259 | A | with pvc's brute-force oracle `brute_segments` (1220) |
| ruzzo_tompa_matches_brute_force_on_pseudo_random | 1278 | A | pvc's 64-bit LCG, `(state >> 33)` as u32; 2000 cases under 160 long |
| finds_a_clean_perfect_tract | 1308 | A | |
| a_substitution_keeps_one_interval | 1324 | A | |
| a_single_indel_keeps_one_interval | 1348 | A | |
| an_n_run_emits_nothing | 1371 | A | |
| a_soft_masked_tract_is_found | 1378 | A | step 1 compares letters regardless of case, as pvc's does |
| only_the_primitive_period_is_emitted | 1388 | A | |
| a_genuine_higher_primitive_period_is_kept | 1405 | A | |
| a_homopolymer_is_only_period_one | 1420 | A | |
| a_period_multiple_re_detection_is_dropped_by_geometry | 1441 | A | |
| the_min_copies_floor_drops_short_segments | 1469 | A | |
| segment_options_default_is_cap_window_no_smoothing | 1196 | C | region scanner's options |
| empty_input_yields_no_regions, a_lone_repeat_tiles_as_unique_repeat_unique, adjacent_repeats_of_different_periods_merge_into_one_region, a_repeat_over_the_cap_is_a_satellite, merge_gap_bridges_two_nearby_repeats, min_repeat_len_reclassifies_a_short_blip_as_unique | 1517–1621 | C | the region scanner `RegionScanner` |
| scan_windowed_* (5), a_window_carries_its_margins_detections_and_says_what_it_fetched, scan_windowed_reports_a_fetch_failure_once_then_stops | 1681–1959 | C | the windowed scan |
| stream_is_window_invariant_* (3), stream_detects_a_satellite_within_a_window | 2072–2099 | C | the windowed stream |

## segment_criteria.rs, 58 tests: 42 ported

| test | line | class | note |
|---|---|---|---|
| a_bundle_member_is_trimmed_to_whole_copies_and_carries_its_motif | 1496 | B | the first half starts the bases at 101, a window; re-based to 0, it expects (2, 26). The second half is A |
| locus_reports_its_span_period_and_purity | 1624 | B | without the contig name |
| is_compound_flags_atat_but_not_at_or_atc | 1665 | A | |
| is_compound_spares_a_unit_dominated_by_one_base_in_every_rotation | 1677 | A | |
| purity_is_one_for_perfect_tiling | 1708 | A | |
| purity_drops_with_one_interruption | 1714 | A | |
| minimal_trim_strips_partial_motifs_at_both_ends | 1719 | A | |
| minimal_trim_returns_none_without_two_clean_copies | 1726 | A | |
| split_bundles_sets_aside_whole_clusters_and_keeps_isolated | 1731 | A | |
| split_bundles_keeps_all_when_none_close | 1750 | A | |
| classify_emits_a_clean_locus_one_based | 1766 | B | params; without the contig name |
| classify_drops_period_over_six_and_compound | 1787 | A | params |
| classify_drops_below_copy_number_floor | 1796 | A | params |
| classify_keeps_a_locus_whose_flank_clamps_to_one_base | 1805 | A | params; step 6's contig-end test, not step 2 |
| each_rejection_reason_names_the_gate_that_fired | 1831 | A | the port keeps the reason of each rejection |
| the_walk_reaches_only_one_of_classifications_five_gates | 1938 | A | |
| classify_drops_locus_with_empty_left_flank | 1977 | A | params |
| classify_drops_locus_with_empty_right_flank | 1987 | A | params |
| classify_drops_period_one_homopolymer_and_spares_the_neighbour_ssr | 1996 | A | params |
| prefilter_drops_low_copy_noise | 2016 | A | params |
| prefilter_keeps_a_higher_period_tract_that_is_not_a_multiple | 2031 | A | params |
| classify_upper_cases_a_soft_masked_tracts_motif | 2318 | B | drop `assert_agrees`, the comparison with pvc's production code |
| classify_reports_the_trimmed_start_not_the_detected_one | 2358 | B | drop `assert_agrees` |
| classify_keeps_an_imperfect_tract_above_the_floor_and_drops_one_below | 2396 | B | drop `assert_agrees_at`; params at purity 0.8, then 0.99 |
| the_copy_floor_counts_the_trimmed_tract_not_the_detected_one | 2452 | A | params, then copy floors (10, 5, 4, 3, 2, 3) |
| classify_applies_the_score_gate_when_the_floor_actually_bites | 2500 | B | drop `assert_agrees_at`; params with score floor 50 |
| classify_at_default_never_rejects_a_score_the_scanner_can_emit | 2559 | B | pvc-default |
| classify_rejects_a_zero_flank | 2930 | B | a bundle distance of 0 raises ValueError, where pvc panics |
| classify_hands_back_the_bundle_members_it_sets_aside | 2952 | B | drop `assert_agrees_at`; `catalog_criteria()` (1525): params with bundle distance 50 |
| classify_bundles_only_what_the_flank_test_rejects | 3007 | A | params |
| narrowing_the_period_scope_changes_what_is_admitted | 3053 | A | params with periods (2,2), then (1,6) |
| classify_bundles_at_the_flank_radius | 3104 | B | drop the two `assert_agrees_at`; bundle distance 50, then 5 |
| the_minimum_is_looked_up_by_the_tracts_own_period | 3159 | A | params with floors (10, 9, 4, 3, 3, 3), then (10, 5, 11, 3, 3, 3) |
| min_copies_uses_the_fallback_only_beyond_the_table | 3216 | B | the port's floors cover periods 1 to 6 only; keep the uniform table of 4 |
| classify_rejects_a_period_ceiling_no_motif_can_hold | 3246 | B | periods (2,7) raise ValueError |
| classify_rejects_a_nan_purity_floor | 3264 | B | a NaN purity raises ValueError |
| raising_the_copy_floor_changes_what_is_admitted | 3284 | A | params, then a floor of 9 at every period |
| joins_cluster_is_strict_at_the_threshold | 3336 | A | |
| a_tract_inside_the_cluster_always_joins_it | 3352 | A | |
| a_member_spanning_past_its_neighbour_keeps_the_cluster_open | 3373 | A | |
| prefilter_drops_an_inverted_interval_instead_of_underflowing | 3404 | A | params |
| differential_vs_production_build_loci_on_scanner_output | 3423 | B | `synthetic_ref.fa`; drop `assert_agrees`; what is left checks that loci are found |
| motif_keeps_bytes_verbatim_and_reports_its_period, motif_rejects_lengths_outside_the_ssr_range, motif_equality_ignores_the_unused_tail | 1565–1592 | C | pvc's `Motif` type; the motif is a `str` here |
| motif_matches_productions_on_the_same_bytes | 1607 | C | pvc's production code |
| locus_rejects_position_zero_because_it_is_one_based, locus_rejects_an_inverted_span, locus_rejects_non_finite_purity | 1636–1656 | C | `SsrSegment::new`'s checks; the port's segment comes only from its own steps |
| differential_vs_production_build_loci_on_crafted_cases | 2207 | C | pvc's production code; its inputs are in the classify tests |
| the_catalog_copy_floor_table_reproduces_both_of_productions | 2616 | C | pvc's production code |
| classify_keeps_a_locus_the_window_edge_would_have_eaten, classify_still_drops_a_locus_at_the_real_contig_end, a_tract_at_the_windows_left_edge_keeps_its_flank_from_the_contig, a_tract_at_the_windows_right_edge_keeps_its_flank_from_the_contig, bundled_coordinates_are_window_offsets_not_genomic | 2697–2843 | C | `classify` over a window; the port always classifies a whole contig |
| classify_rejects_a_zero_bases_start, classify_rejects_a_contig_len_that_understates_the_contig | 2889, 2912 | C | arguments of the windowed `classify` |

## region_typing/mod.rs, 49 tests: 31 ported

| test | line | class | note |
|---|---|---|---|
| absorption_by_sweep_joins_satellites_that_abut_after_growing_at_radius_zero | 1591 | B | keep the fixed point's half: satellites 22–86 and 88–88 (1-based) and a member at 86–88 (0-based) give one satellite 22–88, the member absorbed; drop the sweep's half |
| a_lone_tract_partitions_as_generic_locus_generic | 2037 | A | catalog |
| a_long_array_is_one_satellite_and_swallows_the_locus_inside_it | 2081 | B | pvc-default, then satellite cap 10000 |
| abutting_coverage_runs_merge_into_one | 2140 | A | |
| a_repeat_free_contig_is_one_generic_region | 2214 | B | pvc-default |
| an_impure_tract_is_generic_when_its_pieces_fall_under_the_copy_floor | 2250 | A | catalog |
| a_homopolymer_is_generic_and_does_not_take_its_neighbour_with_it | 2298 | B | catalog; drop the last assertion, a tally |
| a_homopolymer_of_six_or_more_is_a_period_one_locus_at_default | 2359 | B | pvc-default |
| the_mononucleotide_copy_floor_is_exactly_eight | 2441 | B | pvc-default; pvc_eval's own floor of 5 has its own test |
| the_dinucleotide_copy_floor_is_exactly_six | 2465 | B | pvc-default; likewise |
| a_rejected_repeat_leaves_no_hole | 2486 | B | pvc-default |
| a_tract_at_position_one_is_generic_and_the_partition_still_starts_at_one | 2512 | B | pvc-default |
| an_empty_contig_has_no_regions | 2530 | B | pvc-default |
| two_close_tracts_become_one_bundle_carrying_both | 2553 | B | pvc-default |
| a_bundle_longer_than_the_cap_is_a_satellite_and_is_counted_as_one | 2629 | B | pvc-default with each cap; keep the type at the hull, drop the four counts |
| a_bundle_at_a_contig_end_is_dropped_unless_capped | 2679 | B | pvc-default with each cap |
| the_hull_ends_past_a_ragged_last_member | 2798 | B | pvc-default |
| three_chained_tracts_become_one_bundle_of_three | 2847 | A | catalog |
| a_tract_clear_of_a_bundle_is_still_admitted | 2881 | B | pvc-default |
| a_low_copy_repeat_is_generic_not_a_bundle | 2917 | B | pvc-default |
| bundle_clusters_regroups_exactly_what_the_split_set_aside | 2938 | A | |
| the_resident_partition_reproduces_the_golden_catalog | 2977 | B | periods 2 to 6, floors (10, 5, 4, 3, 3, 3), bundle distance 50, purity 0.8, satellite cap 1000, bundle cap 100, scan 2/7/2; both fixtures |
| a_bed_returns_the_same_findings_the_whole_genome_run_does | 3120 | A | catalog; the requested spans are a list |
| a_microsatellite_beside_a_satellite_is_absorbed_into_it | 3230 | B | pvc-default |
| absorption_tracks_the_bundle_radius | 3316 | B | pvc-default with bundle distance 10, 15, 20 and 30 |
| the_satellite_cap_applies_to_the_cleaned_coverage_not_the_raw | 3370 | B | pvc-default |
| a_bed_does_not_change_what_things_are | 3548 | B | pvc-default |
| a_bed_span_is_tiled_by_what_comes_back | 3608 | B | pvc-default |
| generic_clips_at_the_edge_and_findings_come_back_whole | 3639 | B | pvc-default |
| a_locus_straddling_the_bed_edge_is_emitted_whole | 3712 | B | pvc-default |
| two_spans_sharing_a_scan_span_do_not_leak_the_gap_between_them | 3750 | B | pvc-default |
| absorption_by_sweep_is_the_iterations_fixed_point | 1551 | C | property test of the sweep against the fixed point; only the fixed point is ported |
| from_normalized_spans_rebuilds_a_normalized_set_exactly, from_normalized_spans_refuses_what_this_pipeline_never_writes, from_normalized_spans_survives_a_span_ending_at_the_coordinate_ceiling, from_normalized_spans_accepts_every_normalized_set (property test), whole_contigs_covers_each_contig_end_to_end, the_seam_widens_but_does_not_rebase, a_zero_length_contig_is_dropped_before_a_consumer_sees_it, from_bed_path_parses_converts_and_coalesces, a_bad_bed_is_rejected_before_any_region_is_typed | 1636–1939, 3772 | C | pvc's region sets and BED reading |
| overlapping_coverage_counts_a_shared_base_once, counts_start_at_zero, a_bundle_dropped_at_a_contig_end_is_charged_as_flank_clamped, counts_tally_the_regions_yielded | 1776, 1843, 2740, 3416 | C | the tallies |
| default_config_is_the_short_read_settings | 1795 | C | pvc's own defaults |
| a_satellite_merges_with_the_generic_ground_either_side, a_satellite_across_a_bed_edge_is_clipped_to_each_requested_span, ground_without_a_satellite_is_returned_unchanged | 3817–3873 | C | `hand_satellites_to_the_generic_path` |

The tests that take requested spans (`a_bed_*`, `generic_clips_*`, `a_locus_straddling_*`,
`two_spans_sharing_*`) test pvc's `emit_into`. The port has that rule in `str-regions --regions`,
for STR regions, and these tests run against it. Their assertions on generic regions, which pvc cuts
at the edge of a span, are kept only if the port keeps `emit_into` whole as a private step. That is
the plan's choice, and the dropped assertions are listed in the test.

## repeat_catalog/row.rs: 2 ported, for step 2

| test | line | class | note |
|---|---|---|---|
| fifteen_bases_of_flank_is_kept_and_fourteen_is_not | 338 | B | against step 2: (CAG)×6 with 15, then 14, G's before it and 40 after; kept, then dropped |
| the_right_flank_is_measured_against_the_contigs_own_end | 363 | B | the same at the right end |

The other eight tests of row.rs test the catalogue's rows, which the port does not build.
