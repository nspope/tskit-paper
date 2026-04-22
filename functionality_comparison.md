# Functionality comparison: justification

This document is the long-form companion to Table S2 in the paper
(`functionality_table.tex`). For every operation listed in the table,
it records the cell assignment (✓ full, ◐ partial, blank none) for each
of the three comparison libraries — ARGneedle-lib, matUtils (BTE), and
DendroPy — together with justification/references. The
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
  BTE Python interface:
  <https://jmcbroome.github.io/BTE/build/html/index.html>;
  source: <https://github.com/jmcbroome/BTE/blob/main/src/bte.pyx>.

- **DendroPy** — library reference:
  <https://jeetsukumaran.github.io/DendroPy/library/index.html>.

## Library versions and documentation snapshot

The cell assignments below are defended against a frozen snapshot of
each library's user-facing documentation and source, captured on
**2026-04-11** and kept in the working directory `lib_docs/` (not
tracked in git). The pinned versions are:

| Library        | Version | Release date |
|----------------|---------|--------------|
| tskit          | 1.0.2   | (as per `tskit_python-api.html`) |
| arg-needle-lib | 1.2.1   | 2025-11-07 |
| UShER/matUtils | 0.6.6   | 2024-07-11 |
| BTE            | 0.9.3   | (latest tagged release on GitHub at snapshot time) |
| DendroPy       | 5.0.8   | (as per `dendropy_index.html`) |

### Reproducing the snapshot

Each file in `lib_docs/` was saved directly from the upstream source
listed below. To re-create the snapshot, fetch each URL at the
corresponding version tag.

Rendered documentation pages (HTML):

- `tskit_python-api.html` — <https://tskit.dev/tskit/docs/v1.0.2/python-api.html>
- `argneedle_manual.html` — <https://palamaralab.github.io/software/argneedle/manual/>
- `argneedlelib_readthedocs.html` — <https://arg-needle-lib.readthedocs.io/en/latest/>
- `argneedlelib_modules.html` — <https://arg-needle-lib.readthedocs.io/en/latest/modules.html>
- `matutils.html` — <https://usher-wiki.readthedocs.io/en/latest/matUtils.html>
- `dendropy_index.html` — <https://jeetsukumaran.github.io/DendroPy/library/index.html>
- `dendropy_popgenstat.html` — <https://jeetsukumaran.github.io/DendroPy/library/popgenstat.html>
- `dendropy_treecompare.html` — <https://jeetsukumaran.github.io/DendroPy/library/treecompare.html>
- `dendropy_treemodel.html` — <https://jeetsukumaran.github.io/DendroPy/library/treemodel.html>

Raw source (at the version tags listed above):

- `argneedle___init__.py`, `argneedle_convert.py`, `argneedle_grm.py`,
  `argneedle_metrics.py`, `argneedle_pybind.cpp` — from
  <https://github.com/PalamaraLab/arg-needle-lib/tree/v1.2.1> under
  `src/arg_needle_lib/` (Python) and `src/` (C++).
- `bte.pyx`, `bte_index.rst` — from
  <https://github.com/jmcbroome/BTE/tree/v0.9.3> under `src/`
  and `docs/source/` respectively.
- `matutils.rst` — from
  <https://github.com/yatisht/usher/tree/v0.6.6> under `docs/`.

## Inclusion criterion

An operation counts for a library only if it is directly accessible
through that library's Python API. (For matUtils the Python API is
BTE, since matUtils itself is CLI-only.)

- **✓ (full):** documented and accessible, with no substantive
  restriction relative to the tskit equivalent.
- **◐ (partial):** accessible (importable and callable from Python)
  but not clearly documented in the user-facing docs, or documented
  but restricted to a subset of what tskit offers.
- **blank:** not available through the Python API at all.

## Low-level API availability

The cell assignments above score each library's Python API. A separate
question is whether a library also exposes a *documented C or C++
library* that external code can link against directly.

- **arg-needle-lib.** The ReadTheDocs site is titled "arg-needle-lib
  Python API Reference" (`argneedlelib_readthedocs.html`) and its
  toctree contains only Python module pages; there is no C++ API
  reference or header-file documentation.
  The manual notes that C++ code exists but does not document it:

  > "For more advanced users, `arg-needle-lib` contains lower-level
  > C++ functions that can be used by including `arg-needle-lib` as
  > a C++ module. However, we have designed `arg-needle-lib` so that
  > users can in most cases simply use the Python API. Please reach
  > out if you have a particular use case not included in the
  > Python API and would like to discuss options for having it
  > implemented."

- **matUtils / UShER.** The matUtils documentation (`matutils.rst`,
  `matutils.html`) documents CLI subcommands only. The C++ code is
  compiled into the `matUtils` and `usher` binaries; no library-level
  API is documented for external C/C++ consumers.

- **BTE.** BTE is itself a Python extension, not a standalone C++
  library. From `bte_index.rst`:

  > "BTE (Big Tree Explorer) is a Python extension for analysis and
  > traversal of extremely large phylogenetic trees. […] BTE
  > streamlines this process by exposing the heavily optimized MAT
  > library underlying UShER and matUtils to Python."

---

## 1. Data model

This section is about properties of the in-memory data structure
each library exposes — what biological entities can be represented
and how they are accessed in code. File-format concerns are in
section 2 (file formats).

### Lossless ARG representation

Here by "ARG" we mean which segments of genome are inherited by whom:
in other words, the topology and times only.
The mutations/genotypes are covered by the next topics.

- **tskit (✓):** the tree-sequence data model encodes a full
  ancestral recombination graph as a shared-edge structure. See
  [`tskit.TreeSequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence)
  and the data-model overview at
  <https://tskit.dev/tskit/docs/stable/data-model.html>.

- **ARGneedle-lib (✓):** the main data structure of the
  library, `arg_needle_lib.ARG`, stores a full
  recombination graph and supports round-tripping via
  `serialize_arg` / `deserialize_arg`.

- **matUtils/BTE (blank):** the MAT data model is a single
  phylogeny - there is no notion of recombination or alternative
  topologies along a genome.

- **DendroPy (blank):** stores a `Tree` or `TreeList`; no ARG
  data structure.

### Mutations integrated with topology

- **tskit (✓):** the
  [Site](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.SiteTable)
  and
  [Mutation](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.MutationTable)
  tables are members of `TableCollection`; mutations
  are anchored to specific edges and inherited by descendants
  during tree traversal.

- **ARGneedle-lib (✓):** mutations are an intrinsic part of the
  ARG; `generate_mutations`, `get_mutations_matrix`, and
  `get_genotype` operate on mutations attached to ARG edges.

- **matUtils/BTE (✓):** mutation annotation is the entire point of
  the data model - `MATNode.mutations` is a first-class node
  property and the protobuf format encodes per-edge mutations
  natively.

- **DendroPy (◐):** sequence/character data lives in a separate
  `CharacterMatrix` indexed by taxon; mutations are not anchored to
  tree edges (but, they are associated with nodes).

### Integrated reference sequence

- **tskit (✓):**
  [`TreeSequence.reference_sequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.reference_sequence)
  /
  [`has_reference_sequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.has_reference_sequence)
  and
  [`TableCollection.reference_sequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TableCollection.reference_sequence)
  store the reference inside the data structure itself, alongside
  schema metadata.

- **ARGneedle-lib (blank):** no reference-sequence concept.

- **matUtils/BTE (◐):** the reference is supplied as an
  external `--input-fasta` (`-f`) argument to `matUtils summary`,
  `matUtils extract`, and as an argument to `MATree.translate`
  (so, "partial").
  It is not stored inside the protobuf.

- **DendroPy (blank):** no reference concept.

### Individual / ploidy abstraction

- **tskit (✓):** the
  [`IndividualTable`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.IndividualTable)
  groups one or more sample nodes under a single biological
  individual, capturing ploidy directly. This makes diploids and
  family-structured simulations first-class.

- **ARGneedle-lib (◐):** there is no representation of individuals,
  but the Python API has some diploid-aware methods: `exact_arg_grm`
  and `monte_carlo_arg_grm` accept a `diploid=True` flag, and
  `haploid_grm_to_diploid` operates on GRMs; all these work by
  pairing neighbouring sample IDs as a haploid couple. Counted as
  partial because there is no representation in data besides this
  adjacent-IDs convention.

- **matUtils/BTE (blank):** `MATNode` is a single-leaf abstraction;
  no individual or ploidy concept.

- **DendroPy (blank):** `Taxon` and `TaxonNamespace` carry labels
  but there is no individual-vs-genome distinction.

### Pedigree encoding

- **tskit (✓):** the
  [`IndividualTable`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.IndividualTable)
  records `parents` for each individual, so pedigrees are stored
  natively inside the tree sequence.

- **ARGneedle-lib (blank):** no pedigree encoding.

- **matUtils/BTE (blank):** no pedigree encoding.

- **DendroPy (blank):** no pedigree encoding.

### Population abstraction

- **tskit (✓):** the
  [`PopulationTable`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.PopulationTable)
  stores population assignments inside the data structure,
  and it is easy to pull sets of samples associated with each population.

- **ARGneedle-lib (blank):** samples are a flat haploid list; no
  built-in population concept.

- **matUtils/BTE (◐):** matUtils accepts geographic
  region annotations as an external TSV (used by `introduce`),
  so this is a "partial",  but
  the MAT data structure itself has no population concept.

- **DendroPy (◐):** population labels can be encoded as taxon
  labels and passed to `PopulationPairSummaryStatistics` (thus, "partial"),
  but they are not a typed part of the tree data structure.

### Columnar (NumPy) data access

- **tskit (✓):** every column of every Table is a NumPy array
  (`nodes.time`, `edges.left`, `edges.right`, `mutations.node`, ...).

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

- **ARGneedle-lib (blank):** no VCF writer.
  Note that `threads` will output to VCF
  (via `threads vcf`), but AFAICT this needs a `.threads` file,
  and there is no current way to convert an `.argn` file file to `.threads`.

- **matUtils/BTE (✓):** `matUtils extract --write-vcf` and
  `MATree.write_vcf` (`lib_docs/bte.pyx`).

- **DendroPy (blank):** no VCF reader or writer.

### Newick export

- **tskit (✓):** [`Tree.as_newick`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.as_newick)
  and [`TreeSequence.as_newick`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_newick).

- **ARGneedle-lib (✓):** `arg_needle_lib.arg_to_newick`.

- **matUtils/BTE (✓):** `matUtils extract --write-newick`;
  `MATree.get_newick` / `MATree.write_newick`.

- **DendroPy (✓):** `Tree.write(schema="newick")` and
  `Tree.as_string(schema="newick")`.

### FASTA export

- **tskit (✓):** [`TreeSequence.as_fasta`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.as_fasta)
  / [`write_fasta`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.write_fasta)
  and [`TreeSequence.alignments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.alignments).

- **ARGneedle-lib (blank):** no FASTA writer.

- **matUtils/BTE (blank):** BTE uses FASTA for various operations,
  but does not write out to FASTA.

- **DendroPy (✓):** `CharacterMatrix.write(schema="fasta")`.

### NEXUS export

- **tskit (✓):** writer, no reader (note this is "export")

- **ARGneedle-lib (blank):** not supported.

- **matUtils/BTE (blank):** not supported.

- **DendroPy (✓):** first-class NEXUS reader/writer via the
  unified `read` / `write` schemas (NeXML is also supported).

---

## 3. Tree operations

### Iterate trees along a genome

- **tskit (✓):** [`TreeSequence.trees`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.trees)
  iterator and [`TreeSequence.breakpoints`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.breakpoints).

- **ARGneedle-lib (blank):** the C++ side iterates local trees,
  (see well-commented `arg-needle-lib/src/arg_traversal.hpp`),
  but the Python API does not expose a per-tree iterator.

- **matUtils/BTE (blank):** the data model is a single phylogeny,
  not a sequence of trees along a genome; no notion of recombination.

- **DendroPy (blank):** same - single tree or list of trees, no
  genomic positioning.

### Tree traversal (pre-/post-order)

- **tskit (✓):** [`Tree.preorder`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.preorder),
  [`Tree.postorder`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.postorder),
  [`Tree.timeasc`/`timedesc`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.timeasc).

- **ARGneedle-lib (blank):** the C++ `arg_traversal.hpp` machinery
  traverses ARG nodes (`time_efficient_visit`), but no Python-level
  pre/post-order iterator over local trees is exposed.

- **matUtils/BTE (✓):** `MATree.depth_first_expansion` (pre-order)
  and `MATree.breadth_first_expansion`.

- **DendroPy (✓):** `Tree.preorder_node_iter`, `postorder_node_iter`,
  `levelorder_node_iter`, `leaf_node_iter`.

### MRCA and common-ancestor queries

- **tskit (✓):** [`Tree.mrca`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.mrca)
  and [`Tree.tmrca`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.tmrca).

- **ARGneedle-lib (✓):** `arg_needle_lib.most_recent_common_ancestor`;
  also see `arg_needle_lib.kc2_tmrca_sv_stab`
  for ARG-vs-ARG TMRCA comparisons.

- **matUtils/BTE (✓):** `MATree.LCA(node_ids)` (last common ancestor)
  and `MATNode.parent` for traversal upward.

- **DendroPy (✓):** `Tree.mrca(taxa=...)`.

### Branch length and total branch length

- **tskit (✓):** `Tree.branch_length`, [`Tree.total_branch_length`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.total_branch_length).

- **ARGneedle-lib (✓):** per-edge length is `e.parent.height -
  e.child.height` on any `ARGEdge` (for example
  `arg.node(1).parent_edges()[0]`), and `local_volume` /
  `total_volume` give total branch areas across all branches.

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
  but good support generally.

### Tree balance and shape statistics

- **tskit (✓):** [`Tree.colless_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.colless_index),
  [`Tree.sackin_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.sackin_index),
  [`Tree.b1_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.b1_index),
  and [`Tree.b2_index`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.b2_index)
  cover the standard balance/imbalance indices.

- **ARGneedle-lib (blank):** no balance metrics.

- **matUtils/BTE (◐):** `MATree.tree_entropy` reports per-split
  entropy and `count_clades_inclusive` gives clade sizes - usable as
  shape descriptors but not the standard balance indices.

- **DendroPy (✓):** `treemeasure.colless_tree_imbalance`,
  `sackin_index`, `pybus_harvey_gamma`.

### Tree topology enumeration and ranking

- **tskit (✓):** [`Tree.rank`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.rank)
  and [`Tree.unrank`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.unrank)
  map each leaf-labelled tree to a canonical `(shape, label)` integer
  pair and back;
  [`tskit.all_trees`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.all_trees),
  [`all_tree_shapes`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.all_tree_shapes),
  and [`all_tree_labellings`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.all_tree_labellings)
  enumerate all topologies of a given size; and
  [`TopologyCounter`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TopologyCounter)
  tabulates topology frequencies across a tree sequence.

- **ARGneedle-lib (blank):** no such enumeration API.

- **matUtils/BTE (blank):** operates on a single
  mutation-annotated tree; no enumeration of alternative topologies.

- **DendroPy (blank):** no such enumeration API.

---

## 4. ARG/tree editing

### Simplify (sample-restricted history)

- **tskit (✓):** [`TreeSequence.simplify`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.simplify).

- **ARGneedle-lib (blank):** no equivalent.

- **matUtils/BTE (blank):** `MATree` has various
  functions for subtree selection, but `simplify` is an ARG operation,
  and matUtils is not an ARG library.

- **DendroPy (blank):** same.

### Subset by sample or clade

- **tskit (✓):** [`TreeSequence.subset`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.subset)
  and `simplify` with a sample list.

- **ARGneedle-lib (blank):** no Python-level subset operation.

- **matUtils/BTE (✓):** `matUtils extract --samples` /
  `--clade` / `--regex`; `MATree.subtree`, `MATree.get_clade`,
  `MATree.get_regex`, `MATree.get_random`.

- **DendroPy (✓):** `Tree.extract_tree_with_taxa`,
  `extract_tree_without_taxa`, `prune_taxa`,
  `prune_leaves_without_taxa`.

### Union of ARGs

- **tskit (✓):** [`TreeSequence.union`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.union).

- **ARGneedle-lib (blank):** no ARG union operation.

- **matUtils/BTE (blank):** no union of two MATs.

- **DendroPy (blank):** `TreeList` concatenates lists of trees but
  there is no shared-history union.

### Keep / delete genomic intervals

- **tskit (✓):** [`TreeSequence.keep_intervals`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.keep_intervals)
  and [`delete_intervals`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.delete_intervals).

- **ARGneedle-lib (✓):** `arg_needle_lib.trim_arg`

- **matUtils/BTE (blank):** the data model has no notion of
  genomic intervals over which the tree changes.

- **DendroPy (blank):** same.

### Trim flanking regions

- **tskit (✓):** [`TreeSequence.trim`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.trim).

- **ARGneedle-lib (✓):** `arg_needle_lib.trim_arg` provides exactly this operation.

- **matUtils/BTE (blank):** not applicable.

- **DendroPy (blank):** not applicable.

### Extend haplotypes

- **tskit (✓):** [`TreeSequence.extend_haplotypes`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.extend_haplotypes)
  returns a new tree sequence in which the span of each ancestral
  node is extended across adjacent marginal trees wherever the
  relevant parent–child relationship continues to hold, producing
  a more parsimonious edge table without changing the genotypes.
  Introduced in Fritze et al. (2026).

- **ARGneedle-lib (blank):** no similar operation.
  library's ARG editing primitives are restricted to trimming.

- **matUtils/BTE (blank):** not applicable.

- **DendroPy (blank):** not applicable.

---

## 5. Population-genetic statistics

### Nucleotide diversity (π), segregating sites

- **tskit (✓):** [`TreeSequence.diversity`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.diversity)
  and [`segregating_sites`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.segregating_sites);
  branch and site mode, windowed.

- **ARGneedle-lib (blank):** no diversity statistics in the Python API.

- **matUtils/BTE (◐):** `MATree.compute_nucleotide_diversity` returns
  mean pairwise nucleotide differences over the MAT;
  no segregating sites (thus, "partial").

- **DendroPy (✓):** `popgenstat.nucleotide_diversity`,
  `popgenstat.num_segregating_sites`, and
  `popgenstat.average_number_of_pairwise_differences`.

### Tajima's $D$

- **tskit (✓):** [`TreeSequence.Tajimas_D`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.Tajimas_D).

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** no.

- **DendroPy (✓):** `popgenstat.tajimas_d`.

### $F_{ST}$ and divergence

- **tskit (✓):** [`TreeSequence.Fst`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.Fst)
  and [`divergence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.divergence).

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** the popgenstat module computes within-sample
  diversity statistics but does not implement Fst or between-population
  divergence directly.

### $f$-statistics ($f_2$, $f_3$, $f_4$, Patterson's $D$)

- **tskit (✓):** [`f2`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.f2),
  [`f3`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.f3),
  [`f4`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.f4).

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.

### Allele frequency spectrum (one-way and joint)

- **tskit (✓):** [`TreeSequence.allele_frequency_spectrum`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.allele_frequency_spectrum)
  in branch and site mode, single- and multi-population.

- **ARGneedle-lib (blank):**
  This could be done using `bitset_volume_map` or `stab_return_all_bitsets` (branch mode),
  or with `get_mutations_matrix` (site mode),
  but this isn't a first-class method, and with large numbers of samples
  the former would be $O(2^n)$.

- **matUtils/BTE (blank):** no.

- **DendroPy (◐):** `popgenstat.unfolded_site_frequency_spectrum`
  computes the 1D unfolded SFS for a character matrix. Counted as
  partial because there is no joint/multi-population SFS.

### Linkage disequilibrium

- **tskit (✓):** [`TreeSequence.ld_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.ld_matrix)
  and the `LdCalculator` interface.

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.

### Branch-mode statistics

- **tskit (✓):** every statistic in `tskit.TreeSequence` accepts
  `mode="branch"`, computing the corresponding branch-length
  statistic on the trees themselves rather than on observed sites.
  This is the duality property and has no analogue in the comparison
  libraries.

- **ARGneedle-lib (blank):** no,
  although for instance the GRM is computed by putting down mutations
  under the infinite-sites model,
  which approximates branch mode.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.

---

## 6. Ancestry and relatedness

### Inheritance assignment

- **tskit (✓):** [`TreeSequence.link_ancestors`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.link_ancestors)
  (and the underlying [`TableCollection.link_ancestors`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TableCollection.link_ancestors))
  returns an edge table describing, for each sample in a specified
  set, which segments of the genome are inherited from which
  members of a specified set of ancestors. Introduced in Tsambos
  et al. (2023).

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** not applicable — a mutation-annotated
  tree has no notion of multiple local ancestors along a genome.

- **DendroPy (blank):** not applicable.

### IBD segment extraction

- **tskit (✓):** [`TreeSequence.ibd_segments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.ibd_segments).

- **ARGneedle-lib (blank):** the library has no IBD-segment
  extractor in its Python API.

- **matUtils/BTE (blank):** not applicable to a single phylogeny.

- **DendroPy (blank):** not applicable.

### Genetic relatedness matrix

- **tskit (✓):** [`TreeSequence.genetic_relatedness_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness_matrix),
  [`genetic_relatedness`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness),
  [`genetic_relatedness_weighted`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness_weighted).
  [`genetic_relatedness_vector`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness_vector).

- **ARGneedle-lib (✓):** `arg_needle_lib.exact_arg_grm` and
  `monte_carlo_arg_grm` (in `arg_needle_lib.grm`);
  accompanied by `gower_center`, `row_column_center`,
  and `write_grm` for downstream use.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.

### Genetic relatedness matrix-vector product

- **tskit (✓):** [`TreeSequence.genetic_relatedness_vector`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genetic_relatedness_vector)
  computes the product of the branch GRM with one or more vectors
  directly from the tree sequence without materialising the full
  matrix, giving $O(n + n_T \log n)$ scaling rather than $O(n^2)$.

- **ARGneedle-lib (✓):** available as `arg_matmul`
  (documented but not shown by readthedocs).

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.

### Principal components analysis

- **tskit (✓):** [`TreeSequence.pca`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.pca)
  computes a randomised PCA directly from the tree sequence by
  repeatedly applying the branch-mode GRM-vector product, without
  materialising the full relatedness matrix.

- **ARGneedle-lib (blank):** no direct PCA API; downstream PCA
  would require materialising a GRM via `exact_arg_grm` /
  `monte_carlo_arg_grm` and then calling an external linear-algebra library.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.

### Pairwise divergence / coalescence times

- **tskit (✓):** [`TreeSequence.divergence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.divergence),
  [`divergence_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.divergence_matrix),
  and `Tree.tmrca`.

- **ARGneedle-lib (✓):** `arg_needle_lib.distance_matrix` and
  `distance_matrix_v2` give pairwise distances. ARG nodes carry
  times directly so coalescence-time queries are first-class.

- **matUtils/BTE (✓):** in the tree context, basically the same as MRCA, above.

- **DendroPy (✓):** `treemeasure.patristic_distance` gives
  branch-length sums between taxa, and `node_ages` returns
  coalescence times for an ultrametric tree.

### Genealogical nearest neighbours (GNN)

- **tskit (✓):** [`TreeSequence.genealogical_nearest_neighbours`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genealogical_nearest_neighbours).

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.

---

## 7. Mutations and variants

### Variant / genotype iteration

- **tskit (✓):** [`TreeSequence.variants`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.variants)
  and [`genotype_matrix`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.genotype_matrix).

- **ARGneedle-lib (✓):** `arg_needle_lib.get_mutations_matrix` and
  `get_genotype` return mutation/genotype matrices,
  and `arg.mutations()` allows iteration over mutations.

- **matUtils/BTE (✓):** `MATree.get_mutation_samples`,
  `get_mutation`, `count_mutation_types`, and `count_haplotypes`
  enumerate variants and the samples carrying each mutation.

- **DendroPy (✓):** `dendropy.datamodel.charmatrixmodel`
  lets you iterate over characters.

### Haplotype iteration

- **tskit (✓):** [`TreeSequence.haplotypes`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.haplotypes)
  and [`alignments`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.alignments).

- **ARGneedle-lib (blank):** no haplotype iteration is exposed
  (`write_mutations_to_haps` writes a HAPS file but there is no
  Python iterator).

- **matUtils/BTE (✓):** `MATree.get_haplotype` reconstructs the full
  mutation set carried by a sample relative to the reference;
  `count_haplotypes` enumerates unique haplotypes.

- **DendroPy (✓):** `dendropy.datamodel.charmatrixmodel`
  also provides this, effectively.

### Mutation placement / parsimony

- **tskit (✓):** [`Tree.map_mutations`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.map_mutations)
  performs Hartigan parsimony to place mutations on a tree.

- **ARGneedle-lib (✓):** `map_genotype_to_ARG`,
  `map_genotype_to_ARG_diploid`, `mutation_match`, and
  `mutation_best` place genotypes optimally onto an ARG.

- **matUtils/BTE (✓):** parsimony placement is the central operation
  of UShER and matUtils; `MATree.simple_parsimony` and
  `MATree.get_parsimony_score` expose it from BTE.

- **DendroPy (✓):** provides `dendropy.model.parsimony.fitch_up_pass`
  (note: requires bifurcating trees).

---

## 8. Visualisation

### SVG/vector tree drawing

- **tskit (✓):** [`Tree.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.draw_svg)
  and [`TreeSequence.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.draw_svg)
  produce styled, configurable SVG.

- **ARGneedle-lib (blank):** no built-in plotting.

- **matUtils/BTE (◐):** matUtils emits Auspice JSON for
  visualisation in an external viewer; it does not draw SVG itself.

- **DendroPy (✓):** produces TikZ output, which is convertable to SVG.

### ARG visualisation
- **tskit (✓):** [`TreeSequence.draw_svg`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.draw_svg)
  draws all trees along the genome with shared coordinate axes.

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** not applicable.

- **DendroPy (blank):** not applicable.

### ASCII / text tree rendering

- **tskit (✓):** [`Tree.draw_text`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.Tree.draw_text)
  and [`TreeSequence.draw_text`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.draw_text).

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** no.

- **DendroPy (✓):** `Tree.as_ascii_plot` and `Tree.print_plot`;
  `Tree.as_tikz_plot` for vector output.

---

## 9. Metadata and provenance

### Structured metadata with schemas

- **tskit (✓):** [`tskit.MetadataSchema`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.MetadataSchema)
  with JSON, struct, and permissive codecs; every table column has
  its own schema and validation.

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (◐):** `MATNode.annotations` carries clade-level
  annotations and `matUtils annotate` adds named clade labels - an
  unstructured kind of metadata. No user-defined schema.

- **DendroPy (◐):** similar to matUtils, taxa and trees carry free-form
  `annotations` collections, but there is no schema or validation layer.

### Provenance recording and validation

- **tskit (✓):** [`TreeSequence.provenances`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.provenances)
  and the [provenance schema](https://tskit.dev/tskit/docs/stable/provenance.html);
  most operations that produce a new tree sequence appends a provenance record.

- **ARGneedle-lib (blank):** no.

- **matUtils/BTE (blank):** no.

- **DendroPy (blank):** no.
