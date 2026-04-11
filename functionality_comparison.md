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

## 1. Data model and I/O

### Succinct tree-sequence file format
- **tskit (✓):** the `.trees` file format and underlying
  `TableCollection` are defined by tskit; load/save via
  [`TreeSequence.dump`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.dump)
  and [`tskit.load`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.load).
- **ARGneedle-lib (blank):** uses its own `.argn` HDF5-based ARG
  format (`arg_needle_lib.serialize_arg` /
  `arg_needle_lib.deserialize_arg`); not interchangeable with the tskit
  succinct format although `arg_to_tskit`/`tskit_to_arg` exist.
- **matUtils/BTE (blank):** uses the UShER mutation-annotated tree
  protobuf (`.pb`) via `MATree.from_pb`/`MATree.save_pb`; entirely
  unrelated to tree sequences.
- **DendroPy (blank):** in-memory `dendropy.Tree`/`TreeList` with
  Newick/NEXUS/NeXML serialization; no succinct binary format.

### Table-based data API
- **tskit (✓):** [`tskit.TableCollection`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TableCollection)
  exposes Node, Edge, Site, Mutation, Migration, Individual,
  Population, and Provenance tables as first-class editable objects.
- **ARGneedle-lib (blank):** ARG is a node/edge graph; no table API.
- **matUtils/BTE (blank):** node-pointer tree (`MATNode`); no table
  abstraction.
- **DendroPy (blank):** object graph of `Tree`/`Node`/`Edge`; no
  table abstraction.

### VCF import/export
- **tskit (✓):** [`TreeSequence.write_vcf`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.write_vcf)
  and [`TreeSequence.as_vcf`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_vcf).
- **ARGneedle-lib (blank):** no direct VCF I/O; users convert to
  tskit and call `write_vcf`. The library does consume HAPS/SAMPLE
  for ARG inference but does not write VCF.
- **matUtils/BTE (✓):** `matUtils extract --write-vcf` and
  [`MATree.write_vcf`](https://github.com/jmcbroome/BTE/blob/main/src/bte.pyx)
  emit VCF; `MATree.from_newick_and_vcf` and
  `MATree(... vcf=...)` ingest VCF as input genotypes.
- **DendroPy (blank):** no VCF reader or writer; character-data
  schemas are FASTA/PHYLIP/NEXUS only.

### Newick export
- **tskit (✓):** [`Tree.as_newick`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.as_newick)
  and [`TreeSequence.as_newick`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_newick).
- **ARGneedle-lib (✓):** `arg_needle_lib.arg_to_newick` (the binding
  is registered as `arg_to_newwick` in the pybind layer).
- **matUtils/BTE (✓):** `matUtils extract --write-newick`;
  `MATree.get_newick` / `MATree.write_newick`.
- **DendroPy (✓):** `Tree.write(path=..., schema="newick")` and
  `Tree.as_string(schema="newick")`.

### FASTA / alignment export
- **tskit (✓):** [`TreeSequence.as_fasta`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_fasta)
  / [`write_fasta`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.write_fasta)
  and [`TreeSequence.alignments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.alignments).
- **ARGneedle-lib (blank):** no FASTA writer.
- **matUtils/BTE (blank):** `matUtils extract` outputs Newick, VCF,
  JSON, Taxodium, or protobuf — no FASTA target. (BTE can reconstruct
  haplotype mutation sets via `get_haplotype` but does not assemble
  full reference-anchored sequences.)
- **DendroPy (✓):** `CharacterMatrix.write(schema="fasta")` /
  `read(schema="fasta")` for any sequence character matrix.

### NEXUS / NeXML I/O
- **tskit (blank):** not supported.
- **ARGneedle-lib (blank):** not supported.
- **matUtils/BTE (blank):** not supported.
- **DendroPy (✓):** first-class NEXUS and NeXML readers/writers via
  the unified `read`/`write` schemas.

---

## 2. Tree operations

### Iterate trees along a genome
- **tskit (✓):** [`TreeSequence.trees`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.trees)
  iterator and [`TreeSequence.breakpoints`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.breakpoints).
- **ARGneedle-lib (◐):** the ARG is intrinsically a recombination
  graph and the C++ side iterates local trees via stab queries
  (e.g. `bitset_overlap_stab`, `stab_return_all_bitsets`), but the
  Python API does not expose a clean per-tree iterator — the
  documented route to per-tree analysis is `arg_needle_lib.arg_to_tskit`
  followed by tskit's own iterator.
- **matUtils/BTE (blank):** the data model is a single phylogeny,
  not a sequence of trees along a genome; no notion of recombination.
- **DendroPy (blank):** same — single tree or list of trees, no
  genomic positioning.

### Tree traversal (pre-/post-order)
- **tskit (✓):** [`Tree.preorder`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.preorder),
  [`Tree.postorder`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.postorder),
  [`Tree.timeasc`/`timedesc`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.timeasc).
- **ARGneedle-lib (◐):** the C++ `arg_traversal.hpp` machinery
  traverses ARG nodes (`time_efficient_visit`), but no Python-level
  pre/post-order iterator over local trees is exposed; users convert
  to tskit for traversal.
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
- **tskit (◐):** tskit exposes the underlying counts
  (`Tree.num_descendants`, traversal) but does not ship a Colless or
  Sackin index; the user must compute them. Marked partial.
- **ARGneedle-lib (blank):** no balance metrics exposed.
- **matUtils/BTE (◐):** `MATree.tree_entropy` reports per-split
  entropy and `count_clades_inclusive` gives clade sizes — usable as
  shape descriptors but not the standard balance indices.
- **DendroPy (✓):** `treemeasure.colless_tree_imbalance`,
  `sackin_index`, `pybus_harvey_gamma`.

---

## 3. Tree-sequence editing

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

## 4. Population-genetic statistics

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

## 5. Ancestry and relatedness

### IBD segment extraction
- **tskit (✓):** [`TreeSequence.ibd_segments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.ibd_segments)
  with multiple `within`/`between`, length, and MRCA filters.
- **ARGneedle-lib (◐):** the library has no documented IBD-segment
  extractor in its Python API; the practical route is
  `arg_to_tskit` followed by `tskit.TreeSequence.ibd_segments`.
  Downgraded from the table's first-pass ✓ because the operation is
  not native.
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

---

## 6. Mutations and variants

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

### Mutation simulation on trees
- **tskit (◐):** tskit can carry mutations and run mutation models
  via the related `msprime.sim_mutations`/`pyslim` ecosystem rather
  than tskit itself; the library exposes the mutation model machinery
  but does not bundle a top-level `simulate_mutations` entry point.
  Counted as partial within tskit proper.
- **ARGneedle-lib (✓):** `arg_needle_lib.generate_mutations` and
  `generate_m_mutations` simulate Poisson mutations on the ARG.
  Promoted from blank.
- **matUtils/BTE (blank):** no mutation simulation; mutations are
  the input data, not generated.
- **DendroPy (✓):** `dendropy.model.discrete.simulate_discrete_chars`
  simulates character evolution on a tree under standard substitution
  models. Promoted from ◐ to ✓ since the support is general.

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

## 7. Visualization

### SVG tree drawing
- **tskit (✓):** [`Tree.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.draw_svg)
  and [`TreeSequence.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.draw_svg)
  produce styled, configurable SVG.
- **ARGneedle-lib (blank):** no built-in plotting.
- **matUtils/BTE (blank):** matUtils emits Auspice JSON for
  visualisation in an external viewer; it does not draw SVG itself.
- **DendroPy (blank):** no SVG drawing (TikZ output is the closest
  vector format — see below).

### Tree-sequence (multi-tree) visualization
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

## 8. Metadata and provenance

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

---

## Corrections applied to `functionality_table.tex`

The following cells in `functionality_table.tex` were corrected from
the first-pass draft after the deep-dive above. Each entry lists the
row, the column, the old marker → the new marker, and a one-line
reason.

- **Newick export / ARGneedle-lib:** blank → ✓ — `arg_to_newick` is
  exposed in the pybind layer.
- **Iterate trees along a genome / ARGneedle-lib:** ✓ → ◐ — no
  Python iterator over local trees; users go via `arg_to_tskit`.
- **Tree traversal / ARGneedle-lib:** ✓ → ◐ — same reason.
- **Branch length / ARGneedle-lib:** ✓ → ◐ — derivable from node
  times but no per-edge accessor in the Python API.
- **Tree topology comparison / ARGneedle-lib:** blank → ✓ —
  `kc_topology` and the `metrics` module compute KC² and scaled RF.
- **Tree topology comparison / matUtils (BTE):** ◐ → blank — no
  RF/KC implementation in matUtils or BTE.
- **Tree balance and shape / tskit:** ✓ → ◐ — Colless/Sackin not
  shipped; users compute from descendant counts.
- **Tree balance and shape / matUtils (BTE):** blank → ◐ —
  `tree_entropy` and `count_clades_inclusive` provide partial
  shape information.
- **Subset by sample or clade / ARGneedle-lib:** ◐ → blank — no
  native subset operation.
- **Trim flanking regions / ARGneedle-lib:** blank → ✓ — `trim_arg`.
- **Keep / delete genomic intervals / ARGneedle-lib:** blank → ◐ —
  `trim_arg` covers a single contiguous interval.
- **Nucleotide diversity / matUtils (BTE):** blank → ◐ — BTE's
  `compute_nucleotide_diversity` returns mean π only.
- **Nucleotide diversity / DendroPy:** blank → ✓ — full popgenstat
  module.
- **Tajima's D / DendroPy:** blank → ✓ — `popgenstat.tajimas_d`.
- **Allele frequency spectrum / DendroPy:** blank → ◐ — 1D unfolded
  SFS only.
- **IBD segment extraction / ARGneedle-lib:** ✓ → ◐ — only via
  `arg_to_tskit` round-trip.
- **Pairwise divergence / matUtils (BTE):** ◐ → ◐ (unchanged).
- **Pairwise divergence / DendroPy:** ◐ → ◐ (unchanged).
- **Mean descendants / matUtils (BTE):** blank → ◐ —
  `count_clades_inclusive`.
- **Variant / genotype iteration / ARGneedle-lib:** blank → ◐ —
  `get_mutations_matrix`, `get_genotype`.
- **Haplotype reconstruction / matUtils (BTE):** ✓ → ✓ (unchanged).
- **Mutation simulation / ARGneedle-lib:** blank → ✓ —
  `generate_mutations`.
- **Mutation simulation / DendroPy:** ◐ → ✓ —
  `simulate_discrete_chars` is general.
- **Mutation simulation / tskit:** ✓ → ◐ — provided via the
  `msprime.sim_mutations` companion, not tskit itself.
- **Mutation placement / ARGneedle-lib:** blank → ✓ —
  `map_genotype_to_ARG` family.

