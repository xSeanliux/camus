package graphs

import (
	"fmt"
	"slices"
	"strings"

	"github.com/evolbioinfo/gotree/tree"
)

type Network struct {
	NetTree       *tree.Tree        // tree from extended newick
	Reticulations map[string]Branch // reticulation branches
}

// dummy leaf graftHybrid hangs below a hybrid node; Newick() strips it
const placeholderTip = "####"

const (
	Ui = iota // index of u and w in reticulation array
	Wi
)

type Branch struct {
	IDs [2]int // {0: u, 1: w}
}

func (br Branch) Empty() bool {
	return br.IDs == [2]int{0, 0}
}

// Makes extended newick network out of newick tree and branch data computed by
// the CAMUS algorithm. If bidirectional, each branch whose reverse is valid
// (u not an ancestor of w) is written as two hybrids u->w and w->u, both donors
// above both hybrid points.
func MakeNetwork(td *TreeData, branches []Branch, bidirectional bool) *Network {
	td = td.Clone()
	ret := make(map[string]Branch)
	slices.SortFunc(branches, func(br1, br2 Branch) int {
		top1 := td.LCA(br1.IDs[0], br1.IDs[1])
		top2 := td.LCA(br2.IDs[0], br2.IDs[1])
		return td.Depths[top1] - td.Depths[top2]
	})
	h := 0
	for _, branch := range branches {
		u, w := td.IdToNodes[branch.IDs[Ui]], td.IdToNodes[branch.IDs[Wi]]
		h++
		fwd := fmt.Sprintf("#H%d", h)
		ret[fwd] = branch
		if bidirectional && !td.Under(branch.IDs[Ui], branch.IDs[Wi]) {
			h++
			rev := fmt.Sprintf("#H%d", h)
			ret[rev] = Branch{IDs: [2]int{branch.IDs[Wi], branch.IDs[Ui]}}
			graftTip(td, u, fwd)
			graftTip(td, w, rev)
			graftHybrid(td, w, fwd)
			graftHybrid(td, u, rev)
		} else {
			graftTip(td, u, fwd)
			graftHybrid(td, w, fwd)
		}
	}
	cleanTree(&td.Tree)
	return &Network{NetTree: &td.Tree, Reticulations: ret}
}

// grafts tip named name onto n's parent edge; returns the tip
func graftTip(td *TreeData, n *tree.Node, name string) *tree.Node {
	e, err := n.ParentEdge()
	if err != nil {
		panic(fmt.Sprintf("error in MakeNetwork getting parent edge of node %d: %s", n.Id(), err))
	}
	r := td.NewNode()
	r.SetName(name)
	if _, _, _, err := td.GraftTipOnEdge(r, e); err != nil {
		panic(err)
	}
	return r
}

// grafts hybrid node named name onto n's parent edge
func graftHybrid(td *TreeData, n *tree.Node, name string) {
	p, err := graftTip(td, n, placeholderTip).Parent()
	if err != nil {
		panic(fmt.Sprintf("error in MakeNetwork after grafting onto node %d: %s", n.Id(), err))
	}
	p.SetName(name)
}

func (ntw *Network) Newick() string {
	nwk := ntw.NetTree.Newick()
	nwk = strings.ReplaceAll(nwk, placeholderTip+",", "")
	nwk = strings.ReplaceAll(nwk, ","+placeholderTip, "")
	return nwk
}

// Deletes all branch lengths and support values (since they might be misleading)
func cleanTree(tre *tree.Tree) {
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if e != nil {
			e.SetSupport(tree.NIL_SUPPORT)
			e.SetLength(tree.NIL_LENGTH)
		}
		return true
	})
}

func (ntw *Network) Level1(td *TreeData) bool {
	branches := make([]string, 0)
	for k := range ntw.Reticulations {
		branches = append(branches, k)
	}
	for i := range branches {
		for j := i + 1; j < len(branches); j++ {
			r1 := ntw.Reticulations[branches[i]]
			r2 := ntw.Reticulations[branches[j]]
			vR1 := td.LCA(r1.IDs[0], r1.IDs[1])
			vR2 := td.LCA(r2.IDs[0], r2.IDs[1])
			if vR1 == vR2 || illSorted(vR1, vR2, r1, td) || illSorted(vR2, vR1, r2, td) {
				return false
			}
		}
	}
	return true
}

func illSorted(v1, v2 int, r1 Branch, td *TreeData) bool {
	return td.Under(v1, v2) && (td.Under(v2, r1.IDs[0]) || td.Under(v2, r1.IDs[1]))
}
