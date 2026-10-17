# Report: where δ_T(q, (u,w)) lives in CAMUS

Paper: Willson & Warnow, "CAMUS: scalable phylogenetic network estimation" (Bioinformatics 2026, btag245).
Paper context: Definition 4 defines δ_T(q,(u,w)) — the number of times q appears in the multiset Q if q is satisfied by an edge e added above the ordered pair (u,w), else 0. Equation 1 defines s_e(u,w) = Σ_{q∈Q} δ_T(q,(u,w)), precomputed for all pairs before the DP (Section 2.3 says each δ_T evaluation is O(1) after LCA/leafset/depth preprocessing).

## Direct answer

**δ_T(q,(u,w)) is implemented by `QuartetScore` in `internal/score/quartettotals.go:107-157`.**

```go
// Calculates whether a specific quartet is added by a specific edge.
func QuartetScore(q gr.Quartet, u, w, v, wSub *tree.Node, td *gr.TreeData) int
```

It returns `gr.Qeq` when q is satisfied by the edge added above (u,w) (δ_T contributes `td.NumQuartet(q)`, the multiset multiplicity), `gr.Qneq` when q attaches to the cycle but the topology does not match, and `gr.Qdiff` when q is not applicable to this pair (both ⇒ δ_T = 0).

**s_e(u,w) (Eq. 1) is `quartetsTotal` in `internal/score/quartettotals.go:78-93`**, which sums `QuartetScore == Qeq` hits weighted by `td.NumQuartet(q)` (or 1 in `asSet` mode). `CalculateQuartetTotals` (`quartettotals.go:43-61`) is the paper's preprocessing pass that fills the full `quartetTotals[u][w]` table in parallel, gated by `ShouldCalcEdge` (`quartettotals.go:63-65`: w not under u, cycle length > 3, neither endpoint the root).

## How the code maps to the paper's notation

| Paper | Code |
|---|---|
| q | `q gr.Quartet` |
| ordered pair (u, w) | `u, w` (node IDs in `quartetsTotal`; `*tree.Node` in `QuartetScore`) |
| v = lca_T(u, w) | `v := td.LCA(u, w)` (`quartettotals.go:79`) |
| multiset Q with multiplicities | `td.Quartets(v)` (per-vertex bucket) + `td.NumQuartet(q)` |
| cycle γ_(u,w) side containing w | `wSub` from `getWSubtree` (`quartettotals.go:95-104`) |

Mechanics inside `QuartetScore`:

- `uniqueTaxaBelowNodeFromQ(w, q, td)` (`quartettotals.go:160-171`) requires exactly one quartet taxon below w (`Qdiff` otherwise).
- The loop at lines 114-127 computes, for each of q's four taxa, its LCA with the u-side or w-side of the would-be cycle (`cycleNodes`); `dups` (lines 128, 184-194) rejects degenerate placements.
- Lines 131-151 pick, via precomputed depths (`td.Depths`), the taxon whose attachment point is closest to the reticulation, and compare it against `neighborTaxaQ(q, bi)` — the taxon q's own topology pairs with the below-w taxon. Match ⇒ `Qeq` (satisfied), else `Qneq`.

## Preprocessing that makes δ_T O(1)

All built in `MakeTreeData` (`internal/graphs/treedata.go:27-51`):

- LCA table: `calcLCAs` (`treedata.go:130-166`), O(1) lookup via `td.LCA` (`treedata.go:239-241`)
- Leafset bitsets: `calcLeafset` (`treedata.go:107-127`), queried via `td.InLeafset`
- Depths: `calcDepths` (`treedata.go:169-178`)
- Quartets bucketed to the vertices where they can matter: `mapQuartetsToVertices` (`treedata.go:197-222`), read via `td.Quartets(v)` — this is why `quartetsTotal` iterates `td.Quartets(v)` instead of all of Q

## Where s_e feeds the DP (context)

- The scorer interface `Scorer.CalcScore(u,w,td)` (`internal/score/scorers.go:41-45`) wraps the precomputed table: `MaximizeScorer` returns the raw count (paper's s_e); `NormalizedScorer` and `SymDiffScorer` are variants using `CalculateEdgePenalties` (`internal/score/penalty.go`).
- DP core in `internal/infer/main_dp.go`: `DP[v][k]` is the reticulation-restriction matrix M[v,k] (`main_dp.go:22-28`, filled by `solve`, lines 130-166, driven post-order by `RunDP`, lines 90-104). s(v) (Eq. 3) = `scoreAddEdgeK` (lines 178-209); Pairs(v) enumeration via `SubtreePostOrder`/`SubtreePreOrder` (`helpers.go:9-47`). s_dp (Eqs. 4-5) = `scoreEdgesDown` (u = v case, lines 212-246) and `scoreEdgesAcross` (u ≠ v, lines 249-289), with the per-cycle-node off-path M[l(x)]/M[r(x)] terms accumulated incrementally by `cycleDP.update` (lines 40-72) and the k-budget split by `BestSplit`/`FourWayBestSplit` (`helpers.go:51-110`).
- Backtracking: `traceback.go` (`noCycleTrace`/`cycleTrace`/`cycleTraceNode`); `collateResults` (`main_dp.go:106-127`) emits one branch set per reticulation count k.
- Network assembly: `MakeNetwork` (`internal/graphs/net.go:31-68`); output driven from `camus.go`.
- Separate reuse: `ReticulationScore` (`internal/score/score.go:25-72`) applies the same `QuartetScore` logic per gene tree against a finished network (post-hoc support, not part of the DP).

Line numbers verified against working tree on branch `bidirectional-edges` (2026-09-16); `QuartetScore` and `quartetsTotal` read directly, remainder from codebase survey.
