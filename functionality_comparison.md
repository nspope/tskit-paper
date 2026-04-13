# Functionality comparison: justification

This document is the long-form companion to Table S2 in the paper
(`functionality_table.tex`). For every operation listed in the table,
it records the cell assignment (✓ full, ◐ partial, blank none) for each
of the three comparison libraries — ARGneedle-lib, matUtils (BTE), and
DendroPy — together with the actual function/class/CLI subcommand
that supports the assignment and a documentation link. The
comparisons are illustrative, not exhaustive: the goal is to defend
the markers in Table S2 with concrete pointers a reader can verify.

## Library entry points

- **tskit** — Python API reference:
  <https://tskit.dev/tskit/docs/stable/python-api.html>
- **ARGneedle-lib** — manual:
  <https://palamaralab.github.io/software/argneedle/manual/>;
  Python source under
  <https://github.com/PalamaraLab/arg-needle-lib/tree/main/src>
  (the public package re-exports `arg_needle_lib_pybind`,
  `grm`, `metrics`, `serialize_arg`, and `convert`).
- **matUtils (BTE)** — matUtils CLI:
  <https://usher-wiki.readthedocs.io/en/latest/matUtils.html>;
  BTE Python interface source:
  <https://github.com/jmcbroome/BTE/blob/main/src/bte.pyx>.
- **DendroPy** — library reference:
  <https://jeetsukumaran.github.io/DendroPy/library/index.html>.

Markers are deliberately conservative: a partial mark (◐) is used when
the library can perform the operation only on restricted inputs, only
via conversion through another library, or only with substantial
caveats relative to the tskit equivalent.

---

## 1. Data model

This section is about properties of the in-memory data structure
each library exposes — what biological entities can be represented
and how they are accessed in code. File-format concerns are in
section 2 (file formats).

### Lossless ARG representation
- **tskit (✓):** the tree-sequence data model encodes a full
  ancestral recombination graph as a shared-edge structure. See
  [`tskit.TreeSequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence)
  and the data-model overview at
  <https://tskit.dev/tskit/docs/stable/data-model.html>.
- **ARGneedle-lib (✓):** the headline data structure of the
  library; the `arg_needle_lib.ARG` class stores a full
  recombination graph and supports round-tripping via
  `serialize_arg` / `deserialize_arg`.
- **matUtils/BTE (blank):** the MAT data model is a single
  phylogeny — there is no notion of recombination or alternative
  topologies along a genome.
- **DendroPy (blank):** stores a `Tree` or `TreeList`; no ARG
  data structure.

### Mutations integrated with topology
- **tskit (✓):** the
  [Site](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.SiteTable)
  and
  [Mutation](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.MutationTable)
  tables are first-class members of `TableCollection`; mutations
  are anchored to specific edges and inherited by descendants
  during tree traversal.
- **ARGneedle-lib (✓):** mutations are an intrinsic part of the
  ARG; `generate_mutations`, `get_mutations_matrix`, and
  `get_genotype` (in `arg_needle_lib_pybind.cpp`) operate on
  mutations attached to ARG edges.
- **matUtils/BTE (✓):** mutation annotation is the entire point of
  the data model — `MATNode.mutations` is a first-class node
  property and the protobuf format encodes per-edge mutations
  natively.
- **DendroPy (blank):** sequence/character data lives in a separate
  `CharacterMatrix` keyed by taxon; mutations are not anchored to
  tree edges and traversal does not propagate them.

### Integrated reference sequence
- **tskit (✓):**
  [`TreeSequence.reference_sequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.reference_sequence)
  /
  [`has_reference_sequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.has_reference_sequence)
  and
  [`TableCollection.reference_sequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TableCollection.reference_sequence)
  store the reference inside the data structure itself, alongside
  schema metadata.
- **ARGneedle-lib (blank):** no reference-sequence concept anywhere
  in the pybind layer or Python wrappers.
- **matUtils/BTE (blank):** the reference is supplied as an
  *external* `--input-fasta` (`-f`) argument to `matUtils summary`
  and `matUtils extract` — see `lib_docs/matutils.rst` lines 142
  and 246. It is not stored inside the protobuf.
- **DendroPy (blank):** no reference concept.

### Individual / ploidy abstraction
- **tskit (✓):** the
  [`IndividualTable`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.IndividualTable)
  groups one or more sample nodes under a single biological
  individual, capturing ploidy directly. This makes diploids,
  pedigrees, and family-structured simulations first-class.
- **ARGneedle-lib (◐):** there is no individual table, but the
  Python API has diploid-aware helpers — `exact_arg_grm` and
  `monte_carlo_arg_grm` accept a `diploid=True` flag, and
  `haploid_grm_to_diploid` (in `lib_docs/argneedle_grm.py` lines
  26–37) pairs neighbouring sample IDs as a haploid couple.
  Counted as partial because the convention is implicit.
- **matUtils/BTE (blank):** `MATNode` is a single-leaf abstraction;
  no individual or ploidy concept. (The word "individual" appears
  in `lib_docs/bte.pyx` only in unrelated docstrings.)
- **DendroPy (blank):** `Taxon` and `TaxonNamespace` carry labels
  but there is no individual-vs-sample distinction.

### Population structure encoding
- **tskit (✓):** the
  [`PopulationTable`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.PopulationTable)
  stores population assignments inside the data structure, and the
  population-aware statistics (Fst, divergence, joint AFS) consume
  it directly.
- **ARGneedle-lib (blank):** samples are a flat haploid list; no
  population concept in the pybind layer.
- **matUtils/BTE (blank):** matUtils accepts geographic
  region annotations as an external TSV (used by `introduce`) but
  the MAT data structure itself has no population concept.
- **DendroPy (blank):** population labels can be encoded as taxon
  labels and passed externally to `popgenstat`, but they are not
  a typed part of the tree data structure.

### Columnar (NumPy) data access
- **tskit (✓):** every column of every Table is a NumPy array
  (`nodes.time`, `edges.left`, `edges.right`, `mutations.node`,
  …). This enables vectorised bulk operations and underpins much
  of tskit's performance.
- **ARGneedle-lib (blank):** the ARG is a C++ node-pointer object
  exposed through per-node Python wrappers; no NumPy column views.
- **matUtils/BTE (blank):** `MATNode` is a per-node wrapper around
  the underlying C++ object; tree-wide attributes are accessed
  by traversal, not bulk array.
- **DendroPy (blank):** the `Tree`/`Node`/`Edge` graph is a pure
  Python object hierarchy.

---

## 2. File formats

This section is about which named on-disk file formats each library
can read or write, including each library's *native* binary
serialisation format.

### Efficient binary format
- **tskit (✓):** the `.trees` file format
  (kastore-backed) loaded via [`tskit.load`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.load)
  and written via
  [`TreeSequence.dump`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.dump).
- **ARGneedle-lib (✓):** `.argn` HDF5 format via
  `arg_needle_lib.serialize_arg` / `deserialize_arg`.
- **matUtils/BTE (✓):** UShER mutation-annotated tree protobuf
  (`.pb`) via `matUtils extract --write-pb` and
  `MATree.save_pb` / `MATree.from_pb`.
- **DendroPy (blank):** all on-disk formats are text
  (Newick, NEXUS, NeXML, PHYLIP, FASTA).

### VCF export
- **tskit (✓):** [`TreeSequence.write_vcf`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.write_vcf)
  and [`TreeSequence.as_vcf`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_vcf).
- **ARGneedle-lib (blank):** no VCF writer in the pybind or
  wrapper layers; users round-trip via `arg_to_tskit` and call
  `write_vcf` from tskit.
- **matUtils/BTE (✓):** `matUtils extract --write-vcf` and
  `MATree.write_vcf` (`lib_docs/bte.pyx`).
- **DendroPy (blank):** no VCF reader or writer.

### Newick export
- **tskit (✓):** [`Tree.as_newick`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.as_newick)
  and [`TreeSequence.as_newick`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_newick).
- **ARGneedle-lib (✓):** `arg_needle_lib.arg_to_newick` (registered
  in the pybind layer as `arg_to_newwick`).
- **matUtils/BTE (✓):** `matUtils extract --write-newick`;
  `MATree.get_newick` / `MATree.write_newick`.
- **DendroPy (✓):** `Tree.write(schema="newick")` and
  `Tree.as_string(schema="newick")`.

### FASTA export
- **tskit (✓):** [`TreeSequence.as_fasta`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_fasta)
  / [`write_fasta`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.write_fasta)
  and [`TreeSequence.alignments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.alignments).
- **ARGneedle-lib (blank):** no FASTA writer.
- **matUtils/BTE (blank):** the six FASTA mentions in the matUtils
  documentation (`lib_docs/matutils.rst` lines 59, 63, 122, 142,
  186, 246) are *all* references to the input reference FASTA
  consumed by `--translate` and `--write-taxodium`. Line 186
  explicitly tells users to convert a matUtils-emitted VCF to FASTA
  via the external `vcf2fasta` tool. BTE's `bte.pyx` has no FASTA
  writer either.
- **DendroPy (✓):** `CharacterMatrix.write(schema="fasta")`.

### NEXUS export
- **tskit (blank):** not supported.
- **ARGneedle-lib (blank):** not supported.
- **matUtils/BTE (blank):** not supported.
- **DendroPy (✓):** first-class NEXUS reader/writer via the
  unified `read` / `write` schemas (NeXML is also supported).

---

## 3. Tree operations

### Iterate trees along a genome
- **tskit (✓):** [`TreeSequence.trees`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.trees)
  iterator and [`TreeSequence.breakpoints`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.breakpoints).
- **ARGneedle-lib (blank):** the C++ side iterates local trees via
  stab queries (e.g. `bitset_overlap_stab`, `stab_return_all_bitsets`),
  but the Python API does not expose a per-tree iterator. The
  documented route to per-tree analysis is `arg_needle_lib.arg_to_tskit`
  followed by tskit's own iterator, so the operation is not
  natively available.
- **matUtils/BTE (blank):** the data model is a single phylogeny,
  not a sequence of trees along a genome; no notion of recombination.
- **DendroPy (blank):** same — single tree or list of trees, no
  genomic positioning.

### Tree traversal (pre-/post-order)
- **tskit (✓):** [`Tree.preorder`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.preorder),
  [`Tree.postorder`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.postorder),
  [`Tree.timeasc`/`timedesc`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.timeasc).
- **ARGneedle-lib (blank):** the C++ `arg_traversal.hpp` machinery
  traverses ARG nodes (`time_efficient_visit`), but no Python-level
  pre/post-order iterator over local trees is exposed. Users
  convert to tskit for traversal, so the operation is not natively
  available.
- **matUtils/BTE (✓):** `MATree.depth_first_expansion` (pre-order)
  and `MATree.breadth_first_expansion`.
- **DendroPy (✓):** `Tree.preorder_node_iter`, `postorder_node_iter`,
  `levelorder_node_iter`, `leaf_node_iter`.

### MRCA and common-ancestor queries
- **tskit (✓):** [`Tree.mrca`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.mrca)
  and [`Tree.tmrca`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.tmrca).
- **ARGneedle-lib (✓):** `arg_needle_lib.most_recent_common_ancestor`
  in the pybind layer; also `tmrca_mse` for ARG-vs-ARG TMRCA
  comparisons.
- **matUtils/BTE (✓):** `MATree.LCA(node_ids)` (last common ancestor)
  and `MATNode.parent` for traversal upward.
- **DendroPy (✓):** `Tree.mrca(taxa=...)`.

### Branch length and total branch length
- **tskit (✓):** `Tree.branch_length`, [`Tree.total_branch_length`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.total_branch_length).
- **ARGneedle-lib (◐):** branch lengths derive from `ARGNode` time
  attributes and the C++ traversal returns time intervals
  (`local_volume`, `total_volume` give cumulative branch length); no
  per-edge `branch_length` accessor in the Python API.
- **matUtils/BTE (✓):** `MATNode.branch_length` /
  `MATNode.set_branch_length`.
- **DendroPy (✓):** `Edge.length`, `Tree.length()`,
  `Node.distance_from_root`.

### Tree topology comparison (KC, RF)
- **tskit (✓):** [`Tree.kc_distance`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.kc_distance),
  [`TreeSequence.kc_distance`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.kc_distance),
  and [`Tree.rf_distance`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.rf_distance).
- **ARGneedle-lib (✓):** `arg_needle_lib.kc_topology` and
  `arg_needle_lib.metrics.kc2_tmrca_mse_stab` /
  `rf_total_variation_stab` compute KC² and scaled RF between ARGs.
- **matUtils/BTE (blank):** no Robinson-Foulds, KC, or quartet
  distance is exposed in either matUtils or BTE; topology comparison
  is not part of the matUtils workflow.
- **DendroPy (✓):** `treecompare.symmetric_difference` (unweighted RF),
  `weighted_robinson_foulds_distance`, `euclidean_distance`. No KC,
  but RF is the standard metric so this still counts as full support.

### Tree balance and shape statistics
- **tskit (✓):** [`Tree.colless_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.colless_index),
  [`Tree.sackin_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.sackin_index),
  [`Tree.b1_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.b1_index),
  and [`Tree.b2_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.b2_index)
  cover the standard balance/imbalance indices.
- **ARGneedle-lib (blank):** no balance metrics exposed.
- **matUtils/BTE (◐):** `MATree.tree_entropy` reports per-split
  entropy and `count_clades_inclusive` gives clade sizes — usable as
  shape descriptors but not the standard balance indices.
- **DendroPy (✓):** `treemeasure.colless_tree_imbalance`,
  `sackin_index`, `pybus_harvey_gamma`.

---

## 4. ARG/tree editing

### Simplify (sample-restricted history)
- **tskit (✓):** [`TreeSequence.simplify`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.simplify).
- **ARGneedle-lib (blank):** no native equivalent — the ARG retains
  all ancestry. Users convert to tskit and call `simplify` there.
- **matUtils/BTE (blank):** the closest analogue is `extract`-ing a
  subtree, which prunes leaves but does not remove internal nodes
  that no longer contribute history.
- **DendroPy (blank):** `prune_taxa` removes leaves; there is no
  notion of a sample-restricted ancestral graph to simplify.

### Subset by sample or clade
- **tskit (✓):** [`TreeSequence.subset`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.subset)
  and `simplify` with a sample list.
- **ARGneedle-lib (blank):** no Python-level subset operation; users
  convert to tskit. Downgraded from the table's first-pass mark.
- **matUtils/BTE (✓):** `matUtils extract --samples` /
  `--clade` / `--regex`; `MATree.subtree`, `MATree.get_clade`,
  `MATree.get_regex`, `MATree.get_random`.
- **DendroPy (✓):** `Tree.extract_tree_with_taxa`,
  `extract_tree_without_taxa`, `prune_taxa`,
  `prune_leaves_without_taxa`.

### Union of tree sequences
- **tskit (✓):** [`TreeSequence.union`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.union).
- **ARGneedle-lib (blank):** no ARG union operation.
- **matUtils/BTE (blank):** no union of two MATs.
- **DendroPy (blank):** `TreeList` concatenates lists of trees but
  there is no shared-history union.

### Keep / delete genomic intervals
- **tskit (✓):** [`TreeSequence.keep_intervals`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.keep_intervals)
  and [`delete_intervals`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.delete_intervals).
- **ARGneedle-lib (◐):** `arg_needle_lib.trim_arg` restricts an ARG
  to a single contiguous interval, but there is no general
  multi-interval keep/delete.
- **matUtils/BTE (blank):** the data model has no notion of
  genomic intervals over which the tree changes.
- **DendroPy (blank):** same.

### Trim flanking regions
- **tskit (✓):** [`TreeSequence.trim`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.trim).
- **ARGneedle-lib (✓):** `arg_needle_lib.trim_arg` provides exactly
  this operation. Promoted from blank in the first-pass table.
- **matUtils/BTE (blank):** not applicable.
- **DendroPy (blank):** not applicable.

---

## 5. Population-genetic statistics

### Nucleotide diversity (π), segregating sites
- **tskit (✓):** [`TreeSequence.diversity`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.diversity)
  and [`segregating_sites`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.segregating_sites);
  branch and site mode, windowed.
- **ARGneedle-lib (blank):** no diversity statistics in the Python
  API.
- **matUtils/BTE (◐):** `MATree.compute_nucleotide_diversity` returns
  mean pairwise nucleotide differences over the MAT. Counted as
  partial because it is mean π only — no windowing, no branch mode,
  no segregating-sites variant. Promoted from blank in the
  first-pass table.
- **DendroPy (✓):** `popgenstat.nucleotide_diversity`,
  `popgenstat.num_segregating_sites`, and
  `popgenstat.average_number_of_pairwise_differences`. Promoted
  from blank in the first-pass table.

### Tajima's $D$
- **tskit (✓):** [`TreeSequence.Tajimas_D`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.Tajimas_D).
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not computed.
- **DendroPy (✓):** `popgenstat.tajimas_d`. Promoted from blank.

### $F_{ST}$ and divergence
- **tskit (✓):** [`TreeSequence.Fst`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.Fst)
  and [`divergence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.divergence).
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (blank):** the popgenstat module computes within-sample
  diversity statistics but does not implement Fst or between-population
  divergence directly.

### $f$-statistics ($f_2$, $f_3$, $f_4$, Patterson's $D$)
- **tskit (✓):** [`f2`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.f2),
  [`f3`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.f3),
  [`f4`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.f4).
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (blank):** not exposed.

### Allele frequency spectrum (one-way and joint)
- **tskit (✓):** [`TreeSequence.allele_frequency_spectrum`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.allele_frequency_spectrum)
  in branch and site mode, single- and multi-population.
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (◐):** `popgenstat.unfolded_site_frequency_spectrum`
  computes the 1D unfolded SFS for a character matrix. Counted as
  partial because there is no joint/multi-population SFS and no
  branch mode. Promoted from blank.

### Linkage disequilibrium
- **tskit (✓):** [`TreeSequence.ld_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.ld_matrix)
  and the `LdCalculator` interface.
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (blank):** not exposed.

### Branch-mode statistics on trees
- **tskit (✓):** every statistic in `tskit.TreeSequence` accepts
  `mode="branch"`, computing the corresponding branch-length
  statistic on the trees themselves rather than on observed sites.
  This is the duality property and has no analogue in the comparison
  libraries.
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (blank):** not exposed.

---

## 6. Ancestry and relatedness

### Link ancestors
- **tskit (✓):** [`TreeSequence.link_ancestors`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.link_ancestors)
  (and the underlying [`TableCollection.link_ancestors`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TableCollection.link_ancestors))
  returns an edge table describing, for each sample in a specified
  set, which segments of the genome are inherited from which
  members of a specified set of ancestors. Introduced in Tsambos
  et al. (2023).
- **ARGneedle-lib (blank):** no equivalent in the Python API; the
  library's ARG primitives operate on the whole ARG rather than
  restricting to sample-to-ancestor paths.
- **matUtils/BTE (blank):** not applicable — a mutation-annotated
  tree has no notion of multiple local ancestors along a genome.
- **DendroPy (blank):** not applicable.

### IBD segment extraction
- **tskit (✓):** [`TreeSequence.ibd_segments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.ibd_segments)
  with multiple `within`/`between`, length, and MRCA filters.
- **ARGneedle-lib (blank):** the library has no IBD-segment
  extractor in its Python API. The practical route is
  `arg_to_tskit` followed by `tskit.TreeSequence.ibd_segments`,
  but the operation is not natively available.
- **matUtils/BTE (blank):** not applicable to a single phylogeny.
- **DendroPy (blank):** not applicable.

### Genealogical nearest neighbours (GNN)
- **tskit (✓):** [`TreeSequence.genealogical_nearest_neighbours`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genealogical_nearest_neighbours).
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (blank):** not exposed.

### Genetic relatedness matrix
- **tskit (✓):** [`TreeSequence.genetic_relatedness_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness_matrix),
  [`genetic_relatedness`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness),
  [`genetic_relatedness_weighted`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness_weighted),
  [`genetic_relatedness_vector`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness_vector).
- **ARGneedle-lib (✓):** `arg_needle_lib.exact_arg_grm` and
  `monte_carlo_arg_grm` (in `arg_needle_lib.grm`) are headline
  features and are accompanied by `gower_center`, `row_column_center`,
  and `write_grm` for downstream use.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (blank):** not exposed.

### Pairwise divergence / coalescence times
- **tskit (✓):** [`TreeSequence.divergence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.divergence),
  [`divergence_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.divergence_matrix),
  and `Tree.tmrca`.
- **ARGneedle-lib (✓):** `arg_needle_lib.distance_matrix` and
  `distance_matrix_v2` give pairwise distances; `tmrca_mse` and
  `kc_tmrca_vectors` give TMRCA-based comparisons. ARG nodes carry
  times directly so coalescence-time queries are first-class.
- **matUtils/BTE (◐):** distances on a MAT are mutation-count
  parsimony distances along the tree (recoverable via traversal of
  `MATNode.mutations` and branch lengths); no native pairwise
  divergence-time matrix. Counted as partial.
- **DendroPy (◐):** `treemeasure.patristic_distance` gives
  branch-length sums between taxa, and `node_ages` returns
  coalescence times for an ultrametric tree. Counted as partial
  because there is no built-in pairwise distance matrix function
  beyond looping over `patristic_distance` calls.

### Mean descendants
- **tskit (✓):** [`TreeSequence.mean_descendants`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.mean_descendants).
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (◐):** `MATree.count_clades_inclusive` returns the
  number of leaves under each annotated clade — a related but
  coarser per-clade descendant count. Counted as partial.
- **DendroPy (blank):** no direct equivalent (users iterate
  `leaf_iter` per node).

### Extend haplotypes
- **tskit (✓):** [`TreeSequence.extend_haplotypes`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.extend_haplotypes)
  returns a new tree sequence in which the span of each ancestral
  node is extended across adjacent marginal trees wherever the
  relevant parent–child relationship continues to hold, producing
  a more parsimonious edge table without changing the genotypes.
  Introduced in Fritze et al. (2026).
- **ARGneedle-lib (blank):** no equivalent operation; the
  library's ARG editing primitives are restricted to trimming.
- **matUtils/BTE (blank):** not applicable to mutation-annotated
  trees, which lack the multi-tree structure this operation acts
  on.
- **DendroPy (blank):** not applicable.

---

## 7. Mutations and variants

### Variant / genotype iteration
- **tskit (✓):** [`TreeSequence.variants`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.variants)
  and [`genotype_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genotype_matrix).
- **ARGneedle-lib (◐):** `arg_needle_lib.get_mutations_matrix` and
  `get_genotype` return mutation/genotype matrices but the API is
  oriented around mapping genotypes onto an ARG rather than
  iterating per-site Variant objects. Promoted from blank.
- **matUtils/BTE (✓):** `MATree.get_mutation_samples`,
  `get_mutation`, `count_mutation_types`, and `count_haplotypes`
  enumerate variants and the samples carrying each mutation.
- **DendroPy (blank):** character matrices store sequences but there
  is no per-variant iterator.

### Haplotype reconstruction
- **tskit (✓):** [`TreeSequence.haplotypes`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.haplotypes)
  and [`alignments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.alignments).
- **ARGneedle-lib (blank):** no haplotype iteration is exposed
  (`write_mutations_to_haps` writes a HAPS file but there is no
  Python iterator).
- **matUtils/BTE (✓):** `MATree.get_haplotype` reconstructs the full
  mutation set carried by a sample relative to the reference;
  `count_haplotypes` enumerates unique haplotypes.
- **DendroPy (blank):** sequences are stored verbatim in
  `CharacterMatrix`; no reconstruction from a tree.

### Mutation placement / parsimony
- **tskit (✓):** [`Tree.map_mutations`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.map_mutations)
  performs Hartigan parsimony to place mutations on a tree.
- **ARGneedle-lib (✓):** `map_genotype_to_ARG`,
  `map_genotype_to_ARG_diploid`, `mutation_match`, and
  `mutation_best` place genotypes optimally onto an ARG. Promoted
  from blank.
- **matUtils/BTE (✓):** parsimony placement is the central operation
  of UShER and matUtils; `MATree.simple_parsimony` and
  `MATree.get_parsimony_score` expose it from BTE.
- **DendroPy (blank):** no parsimony placement (DendroPy targets
  Bayesian/likelihood workflows externally).

---

## 8. Visualisation

### SVG tree drawing
- **tskit (✓):** [`Tree.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.draw_svg)
  and [`TreeSequence.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.draw_svg)
  produce styled, configurable SVG.
- **ARGneedle-lib (blank):** no built-in plotting.
- **matUtils/BTE (blank):** matUtils emits Auspice JSON for
  visualisation in an external viewer; it does not draw SVG itself.
- **DendroPy (blank):** no SVG drawing (TikZ output is the closest
  vector format — see below).

### Tree-sequence (multi-tree) visualisation
- **tskit (✓):** [`TreeSequence.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.draw_svg)
  draws all trees along the genome with shared coordinate axes.
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not applicable.
- **DendroPy (blank):** not applicable.

### ASCII / text tree rendering
- **tskit (✓):** [`Tree.draw_text`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.draw_text)
  and [`TreeSequence.draw_text`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.draw_text).
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (✓):** `Tree.as_ascii_plot` and `Tree.print_plot`;
  `Tree.as_tikz_plot` for vector output.

---

## 9. Metadata and provenance

### Structured metadata with schemas
- **tskit (✓):** [`tskit.MetadataSchema`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.MetadataSchema)
  with JSON, struct, and permissive codecs; every table column has
  its own schema and validation.
- **ARGneedle-lib (blank):** ARG nodes carry numeric attributes
  (time, IDs) but no general metadata-schema mechanism.
- **matUtils/BTE (◐):** `MATNode.annotations` carries clade-level
  annotations and `matUtils annotate` adds named clade labels — a
  fixed, unstructured kind of metadata. No user-defined schema.
- **DendroPy (blank):** taxa and trees carry free-form
  `annotations` collections, but there is no schema or validation
  layer.

### Provenance recording and validation
- **tskit (✓):** [`TreeSequence.provenances`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.provenances)
  and the [provenance schema](https://tskit.dev/tskit/docs/stable/provenance.html);
  every operation that produces a new tree sequence appends a
  validated provenance record.
- **ARGneedle-lib (blank):** not exposed.
- **matUtils/BTE (blank):** not exposed.
- **DendroPy (blank):** not exposed.
