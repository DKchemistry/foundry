---
orphan: true
---

# Design: N-to-C macrocycle generation in RFD3

**Status:** Frozen for implementation after user agreement and selected-chain scope review.
**Branch:** `feat/rfd3-cyclic-positional-encoding`
**Reviewed baseline:** `b02eed6a6bdf8f44d14a80cc36e3da13c9f2291c`
**Proof of concept:** `explore/rfd3-binder-trace` at
`7ee29188f90af7c327df1be444e0309a495e45ce`, including
`models/rfd3/docs/examples/macrocycle_binder_poc.md`.

The user has agreed to the public interface, internal transport, and scope below.
Validation is local to the selected cyclic chain; unrelated components retain
existing RFD3 behavior. Implementation and agent review are now authorized.

## 1. Objective

Expose cyclic residue relative positional encoding through RFD3's normal inference
input specification. Support de novo macrocycle monomers and macrocycle binders
using existing checkpoints. Make the change small, explicit, testable, and suitable
for an upstream PR.

The immediate workflows are fixed-length 10- and 12-residue monomer generation and
peptide binder generation against protein targets. Scientific evaluation belongs
to the user. Engineering acceptance concerns input handling, feature correctness,
compatibility, reproducibility of settings, and documentation.

[RFpeptides](https://doi.org/10.1038/s41589-025-01929-w)
([accessible full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12643943/)) describes
cyclic positional encoding for RF2 and RFdiffusion, including monomer and binder
generation. For binders, it changes intrapeptide offsets while retaining the target
and interchain encodings. Its monomer study includes 10- and 12-residue backbones.
This motivates the interface and pair masking here; it does not establish RFD3's
scientific performance.

The local proof of concept wraps residue offsets for hard-coded `asym_id == 0`.
Its report describes visual closure in eight 12-residue insulin receptor binder
outputs. This draft uses that implementation as the behavioral starting point.

## 2. Repository requirements

[CONTRIBUTING.md](../../../CONTRIBUTING.md) requires meaningful names, Google-style
docstrings, PEP 8/20 conventions, tests, logical commits, conventional commit
messages, and early draft PR feedback. It recommends a PR below roughly 400 lines.
Preserve the dependency direction: Foundry depends on AtomWorks for structure
processing; AtomWorks must not acquire an RFD3 dependency.

The guide's statement that tests are unsupported is stale relative to this
checkout. `pyproject.toml`, `.github/workflows/test.yaml`, and
`models/rfd3/tests/conftest.py` define CPU tests and mypy checks. Older cluster-only
tests are excluded explicitly. Add portable tests under `models/rfd3/tests/`.
Ruff uses 88 columns and import sorting. Model documentation is connected through
`models/rfd3/docs/index.rst` and the Sphinx documentation tree.

Keep model-specific policy under `models/rfd3`. No new dependency, checkpoint,
model parameter, or general graph-topology framework is expected. Reuse AtomWorks
annotations and token utilities through the existing RFD3 adapters.

At draft creation, five implementation files already contain uncommitted tracing
changes, and `.trace.md` is untracked. Preserve these during planning. Account for
them explicitly before preparing feature commits; do not stage them accidentally.
The existing remotes are `origin` = `DKchemistry/foundry` and `upstream` =
`RosettaCommons/foundry`. A new fork is not currently needed.

## 3. Agreed first-release scope

Support **one complete de novo canonical amino-acid peptide chain** per design,
using N-to-C cyclic residue relative positional encoding in dialect 2:

- A peptide monomer, with fixed or sampled length. Provide 10- and 12-residue examples.
- A peptide binder against protein targets, using the existing target and hotspot
  conditioning. The peptide may appear before or after the target in the contig.

Motif-containing rings, cyclization of atomized/noncanonical polymers, nucleic
acids or ligands, partial diffusion, symmetry-specific workflows, and multiple
independent macrocycles are future work. These composition restrictions apply to
the selected cyclic chain. Unrelated components elsewhere in an otherwise valid
design retain existing RFD3 behavior, including cofactors and unindexed motifs.
Keep validation proportional; do not build a general topology validator. The
partial-diffusion and active-symmetry exclusions remain request-level mode checks.

Remove the previous proposed minimum of three residues. The offset computation has
no such requirement: it is defined for any nonempty selected chain. Empty or missing
chains are invalid selections. Scientific suitability of a length is not an input
validation rule. Lengths one and two can exercise the arithmetic without implying
that they are useful macrocycle designs.

Training support, sequence design changes, relaxation, closure scoring, benchmark
orchestration, and downstream prediction changes are outside this PR.

### Agreed topology and output behavior

Implement cyclic residue positional encoding. Preserve existing `token_bonds`,
input bonds, terminal features, and structure-writing behavior. Do not add a
special terminal C–N bond to model inputs or written CIF files. Preserve the user's
`cyclic_chains` request in the normal serialized specification/output metadata.

Document that the setting requests cyclic encoding. Geometric closure and
scientific quality are evaluated separately by the user. The flag is not a closure
success measurement. The explicit bond question is resolved for this PR.

## 4. Public interface

Use RF3's public field name on `DesignInputSpecification`:

```python
cyclic_chains: list[str] | None = Field(None, ...)
```

The list contains **assembled design chain IDs**, is case-sensitive, and selects
whole chains. It does not contain input-file chain IDs, integer `asym_id` values,
residue selections, or symmetry names. Omission, null, and `[]` mean linear
encoding, following RF3's optional-list convention. Reject a bare string,
non-string elements, empty IDs, and duplicate IDs with
clear validation errors. Accept at most one selected chain for this release; keep
the list interface consistent with RF3. Do not add a second enable/disable boolean.

Assembled chains currently start at `A`; `/0` advances to the next chain. In
`12,/0,E6-155`, the generated peptide is `A` and the source target segment becomes
`B`. Hotspot selections still refer to the source structure, such as `E64`.
In `E6-155,/0,12`, the generated peptide is `B` and must be selected as `B`.
Document this distinction next to the field and binder example.

Proposed monomer input:

```json
{
  "macrocycle_10": {"length": "10", "cyclic_chains": ["A"]},
  "macrocycle_12": {"length": "12", "cyclic_chains": ["A"]}
}
```

Proposed binder input, using the existing example asset:

```json
{
  "insulinr_macrocycle_12": {
    "dialect": 2,
    "input": "../input_pdbs/4zxb_cropped.pdb",
    "contig": "12,/0,E6-155",
    "cyclic_chains": ["A"],
    "infer_ori_strategy": "hotspots",
    "select_hotspots": {"E64": "CD2,CZ"}
  }
}
```

Paths above assume an example JSON under `docs/examples/`. Publish separate monomer
and binder example files so monomer use does not load target structures. Keep the
usual sampler, seed, checkpoint, and batch controls; this feature needs no special
CLI or Hydra switch. These proposed inputs are not executable until implemented.

## 5. RF3 precedent and transport decision

### 5.1 What RF3 actually does

RF3 is the primary in-repository precedent. Its current implementation is more than
an ID-only pass-through:

1. [`InferenceInput.cyclic_chains`](../../rf3/src/rf3/utils/inference.py) stores the
   public request. `to_pipeline_input()` calls `cyclize_atom_array()`, which moves
   the last carbon near the first nitrogen when needed and adds terminal bonds.
2. [`AddCyclicBonds`](../../rf3/src/rf3/data/cyclic_transform.py) detects terminal
   proximity, updates `token_bonds`, and writes a list of `cyclic_asym_ids` to `feats`.
   Thus RF3 can also recognize cyclic geometry without an explicit chain request.
3. [`RelativePositionEncoding`](../../rf3/src/rf3/model/layers/pairformer_layers.py)
   reads `cyclic_asym_ids`, derives each chain length by counting unique residue
   indices, and chooses the shortest of `d`, `d + L`, and `d - L`. It preserves the
   original offset on equal-distance ties and applies wrapping only within the
   selected chain. Its sign convention is the same `i-j` convention used in RFD3.

Reuse the public name, internal ID name, selection by asym identity, and signed
wrapping semantics. RF3's coordinate/bond-producing steps cannot be reused directly:
RFD3 initializes generated coordinates at the origin, and this PR explicitly excludes
bond conditioning and bond serialization. Resolve the explicit request to IDs at
featurization instead. Do not infer cyclic intent from generated coordinates.

### 5.2 Comparison of representations

| Representation | Benefits | Costs and fit for this scope |
| --- | --- | --- |
| `cyclic_asym_ids`, length derived in RPE | Matches RF3's model-facing vocabulary; compact; uses the existing specification and asym mapping; no new atom annotation. | Needs a small chain-membership reduction in the RPE. For one canonical chain this is straightforward tensor code. Recommended. |
| Per-atom `cyclic_chain_length`, then token feature | Carries a precomputed length and can help detect chain truncation; straightforward broadcasting in RPE. | Adds annotation creation, preservation, padding, defaults, and consistency obligations. The length duplicates information already present in valid canonical token features. No inspected RFD3 requirement justifies it here. |
| Per-token cyclic marker, length derived in RPE | Naturally follows token ordering and uses tensor operations; could be useful for future token filtering. | Introduces another internal convention while still requiring chain resolution and length derivation. A marker is easily derived from RF3-style IDs when needed. No present advantage warrants divergence. |

**Agreed: use `cyclic_chains -> cyclic_asym_ids -> tensor-derived length in RPE`.**
Use RFD3's assembled chain namespace and resolve the IDs where final token features
are available. Retire the proposed `cyclic_chain_length` annotation.

### 5.3 Concrete RFD3 findings

- `DesignInputSpecification.to_pipeline_input()` already places the serialized
  request in `data["specification"]`. That dictionary is available during token
  encoding; `SubsetToKeys` runs at the end of the pipeline.
- `EncodeAF3TokenLevelFeatures` has token-level `chain_id`, `pn_unit_iid`, and the
  computed `asym_id` in the same method. It can map the requested assembled chain
  directly to the actual model identity. Alphabetical position is not a safe proxy.
- The existing pipeline converts `feats` with `ConvertToTorch`, retains the feature
  dictionary through aggregation, and passes it to both RPE instances in
  `TokenInitializer`. No new atom annotation or engine parameter is needed.
- RFD3's `strip_f()` inspects `.shape` on every feature for classifier-free guidance.
  Supply `cyclic_asym_ids` as an integer NumPy array before conversion and an integer
  tensor at the model boundary, rather than relying on transport of a Python list.
  This preserves RF3's meaning and name while satisfying RFD3's feature convention.
- `RFD3.forward()` calls the token initializer before the diffusion rollout, and
  again for the reference features when guidance is enabled. The current optional
  compilation targets are diffusion submodules, not the token initializer. There
  is no evidence that RPE execution requires a stored length for compilation or
  repeated denoising performance. Deriving length with a tensor sum also avoids
  copying RF3's dynamic `unique()` operation.

These are source-level findings. The shell Python inspected during this review has
no Torch, AtomWorks, or Biotite installed, so runtime transport is not yet verified.
The implementation must include a real pipeline and guidance test.

## 6. Implementation contract

### 6.1 Public validation and assembled-chain resolution

Add the field and ordinary schema checks. Reject nonempty cyclic requests with
partial diffusion or active symmetry, and reject multiple selected chains. Because
legacy dialect 1 accepts extras, reject nonempty cyclic requests in both dialect
entry points (`safe_init()` and `create_atom_array_from_design_specification()`).

At assembly, validate that the selected chain exists and is entirely newly generated
canonical peptide. Use `src_component` provenance and existing annotations so a
source-derived motif cannot pass merely because its sequence/coordinates were
unfixed. Apply composition, atomization, unindexing, and continuity checks only
to the selected chain. Do not reject a request because an unrelated retained
component is a ligand/cofactor, nucleic acid, unindexed motif, or noncanonical/
atomized polymer. Existing RFD3 input restrictions still apply to those components.

Require the selected chain to map exclusively to one model asym identity: every
token with that `asym_id` must belong to the selected chain. This prevents a shared
PN unit from including another component in the length reduction or cyclic pair
mask. Reject an ambiguous shared identity instead of guessing a subset or modifying
AtomWorks PN-unit assignment.

Keep these checks in a small helper or the existing validation flow. Avoid custom
bond-graph traversal, a generic topology class, or arbitrary size thresholds.

### 6.2 Feature transport

In `EncodeAF3TokenLevelFeatures`, read the request from
`data.get("specification", {}).get("cyclic_chains") or []`. Ordinary training data and
older inputs lack the request and remain linear. Resolve the selected assembled
`chain_id` using token starts and the encoder's actual `asym_id` array.

For an active request, require a nonempty selected token set, exactly one asym
identity belonging exclusively to that chain, canonical nonatomized residue tokens,
and consecutive distinct residue
indices with one token per residue. These are the invariants needed for the length
calculation, not a general structural validation system. Run checks on residue
identity before the encoder masks sequence names to `GAP`. Target residue numbering
and tokenization need not be contiguous merely because the peptide is cyclic.

Create `feats["cyclic_asym_ids"]` as an integer NumPy array of shape `[1]` for an
active request. For inactive inputs, prefer omitting the feature to leave their
feature dictionaries unchanged; the RPE also accepts an empty integer tensor.
After existing conversion, an active request is an integer tensor of shape `[1]`
on the feature device. Do not expose additional public IDs or lengths.

No constant/annotation registry change is needed. No pipeline reordering is planned.
Verify that the existing parser and padding preserve the normal chain identities
used in resolution. Do not add a second copy of chain lengths solely to audit
hypothetical future cropping: the standard inference route skips training crops.
A real pipeline test must establish that the requested peptide reaches the encoder
in full, and a missing selected chain must raise an error.

Guidance retains cyclic identity in both conditional and reference features; it is
not a feature to zero through `cfg_features`. Unindexed components elsewhere are
appended after indexed tokens by `UnindexFlaggedTokens`; guidance crops that suffix
without renumbering retained asym IDs. A valid selected peptide is indexed and
remains intact. Test the actual `strip_f()` path with an unrelated unindexed suffix,
as well as the one-token edge case because its shape equals the ID-vector length.

### 6.3 RPE arithmetic

Modify `RelativePositionEncodingWithIndexRemoval` locally; do not import RF3's model
class. RFD3's extra unknown-index bin and unindexing override must remain intact.
A small private offset helper is appropriate. Sharing an abstraction across models
is unnecessary for this PR.

For the one selected ID, derive `selected = (asym_id == cyclic_asym_ids[0])` and
`L = selected.sum()`. Keep `L` as an integer tensor. Counting tokens equals counting
residues under the checked canonical, one-token-per-residue contract. This is the
specific reason to use a sum rather than RF3's unique-residue counting. If broader
tokenization support is added later, this assumption must be revisited.

For `d = residue_index[i] - residue_index[j]`, use the shortest signed offset:

```text
wrapped(d, L) = d - L   if 2*d > L
                d + L   if 2*d < -L
                d       otherwise
```

For a valid contiguous chain this is equivalent to RF3's three-candidate selection,
including its tie handling. Apply the result only to `selected[i] & selected[j]`.
Do not use integer ID zero as a special value: zero is a valid cyclic asym ID.
The public feature contract contains at most one ID, so this does not implement
multiple independent macrocycles by treating them as one chain.

Use tensor operations without `.item()`, device-to-host copies, or tensor-value
Python branches. Branching on feature presence or its zero/one element shape is
sufficient. For absent or empty cyclic IDs, take the existing linear calculation.
Do not mutate input tensors. Preserve the existing sequence:

```text
signed residue offset -> optional wrap -> shift/clip/interchain sentinel
                      -> unindexing override -> one-hot -> learned projection
```

Preserve token offsets, entity/chain features, bin counts, state-dict keys and
shapes, and strict checkpoint loading. Both RPE instances receive the same feature.

**Even-length ties:** Keep the ordinary sign at `d = +/-L/2`, matching RF3 and the
proof of concept. For lengths 10 and 12 the ties are respectively `+/-5` and `+/-6`.
For `L > 2`, first-to-last becomes `+1` under RFD3's `i-j` convention and the reverse
becomes `-1`. At length two those pairs are themselves ties; retain their original
sign. At length one the diagonal is zero. These edge cases require no minimum
length rule.

This convention preserves antisymmetry. Opposite pairs in even rings depend on the
chosen sequence origin, so do not claim complete cyclic-relabeling invariance at
those ties. Avoid modulo wrapping that maps both half-ring directions to one sign.

### 6.4 Serialization

Use `get_dict_to_save()` and the ordinary output path to retain non-default
`cyclic_chains`. Test output JSON and `from_rfd3_out()` reload. Keep IDs and derived
lengths internal. Existing sampled-contig metadata records sampled composition;
reloading a range alone is not a promise of resampling the same structure.

No special CIF connectivity, new closure flag, or structure-writer feature is
required. A normal output smoke test is enough to check that the option coexists
with existing writing behavior.

## 7. Verification and acceptance criteria

Use three complementary levels: a small offset-helper test, actual RPE forward
tests, and real input-to-feature-to-RPE integration. Do not expose production
preprojection tensors or add a debug API solely for tests. An ordinary test hook
on the existing linear layer, or independently constructed expected inputs to that
layer, can inspect bins when necessary without changing the production interface.

| Area | Required evidence |
| --- | --- |
| Defaults | Omitted/null/empty public selection and absent/empty cyclic IDs preserve the existing linear RPE output exactly. Input feature tensors are not mutated. |
| Arithmetic | Hand-verified signed offsets for odd and even lengths, including 10/12; diagonal, seam adjacency, shortest distances, antisymmetry, and exact ties. Tiny cases confirm no arbitrary minimum; a larger ring exercises clipping. |
| Pair isolation | Only the selected intrachain residue block changes. Cover binder-first and target-first contigs and a nonzero selected asym ID. |
| Existing encodings | Token offsets, entity/chain terms, interchain sentinel, and unindexing precedence remain unchanged. Use a synthetic unindexing case at model level without supporting motif rings through the public API. |
| Scope checks | Invalid type/ID, missing chain, multiple selections, selected source-derived or atomized residues, and unsupported modes or selected-chain composition fail clearly. Unrelated components are accepted; a shared asym identity is rejected. Reuse fixtures and parameterization; no exhaustive graph-topology suite. |
| Transport | Public input -> assembly -> parsing/padding -> token identity resolution -> integer ID tensor -> actual RPE. Cover monomers and a small protein binder fixture in both contig orders; verify resolved sampled length. Include a positive binder-with-cofactor case and unchanged context/interchain features. |
| Guidance | Actual `strip_f()` preserves the ID tensor and peptide block after removal of an unrelated unindexed suffix, including the one-token shape edge case. |
| Compatibility | State-dict names/shapes remain identical and strict loading succeeds. Existing token/bond features are preserved. |
| Serialization | The public request survives normal output JSON and reload; default omission is preserved. |

Drop multiple-ring support tests, minimum-length rejection tests, annotation-survival
checks for the discarded length field, and generic truncation infrastructure.

Run focused tests, then repository CPU tests, mypy, Ruff checks, and a Sphinx build
for the documentation. Report unavailable dependencies and pre-existing failures
separately. Format changed files only during development to preserve other work.

Later, run a checkpoint/GPU smoke test for a 10-mer, 12-mer, and protein binder,
including a diffusion batch larger than one and the optional compiled path where
available. Verify completion, finite coordinates, chain lengths, and saved settings.
Record the environment/checkpoint. Scientific closure rates and quality remain
the user's evaluation and are not unit-test thresholds.

## 8. Expected change surface and delivery

Paths in this table are relative to `models/rfd3/`.

| File | Responsibility |
| --- | --- |
| `src/rfd3/inference/input_parsing.py` | Public field, compact scope checks, legacy rejection. |
| `src/rfd3/transforms/util_transforms.py` | Resolve assembled chain to cyclic asym ID at token encoding. |
| `src/rfd3/model/layers/blocks.py` | Derive length and conditionally wrap residue offsets. |
| `tests/test_rfd3_cyclic_encoding.py` | Arithmetic, RPE isolation, compatibility, and guidance tests. |
| `tests/test_rfd3_cyclic_inputs.py` | Input, pipeline integration, and serialization tests. |
| `docs/input.md`, `docs/index.rst`, `docs/examples/` | Field reference and runnable monomer/binder examples. |

No change is planned to RF3, `constants.py`, pipeline construction, engine plumbing,
shared Foundry code, sampling, or checkpoint configuration. Revisit only if an
actual integration failure demonstrates a need.

With the transport and narrow validation contract agreed, delegate
input/feature transport and encoder/tests with separate file ownership. An
independent reviewer checks the contract; the coordinating agent owns integration,
documentation consistency, and final review.

Prepare logical conventional commits. Keep existing tracing changes separate and
preserved. Publish the feature branch to the existing personal fork for scientific
tests and record the tested commit. Prepare an upstream draft PR against the
maintainers' current base with the actual validation results.

Aim for the contribution guide's roughly 400-line recommendation without sacrificing
readability or useful tests. This working design can remain personal-branch material;
upstream user documentation should explain the shipped feature concisely.

## 9. Frozen decision record

- Follow RF3 naming: public `cyclic_chains`, model-facing `cyclic_asym_ids`.
- Resolve the explicit request at token encoding and derive length with a tensor
  sum in the RPE; do not add `cyclic_chain_length` annotations.
- Support one complete de novo canonical peptide monomer or protein binder.
  Validate the selected chain's composition and representation, not unrelated
  context. Require its asym identity to exclude every other component.
- Keep partial diffusion, active symmetry, and multiple cyclic selections outside
  this release. Preserve existing context behavior, including cofactors and
  unindexed motifs elsewhere.
- Preserve signed shortest offsets and original-sign half-ring ties as in RF3 and
  the proof of concept. No arbitrary minimum ring length.
- Keep explicit input/output terminal bonds and changes to `token_bonds` outside
  this PR. Serialize the public request normally.
- Use helper, actual RPE, pipeline, guidance, and serialization tests without new
  production debugging APIs. Scientific quality remains the user's evaluation.

The scope review found no conflict requiring a blanket restriction on unrelated
components. The user authorized freezing the design and proceeding in that case.
There are no outstanding design questions. Runtime validation remains implementation
work; record any concrete failure and fix it within this contract where possible.
