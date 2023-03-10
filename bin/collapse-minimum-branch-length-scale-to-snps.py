#!/usr/bin/env python

import argparse
import treeswift

def main(args):
    
    t = treeswift.read_tree(args.tree, 'newick')

    t.collapse_short_branches(args.min_branch_length)

    # Leaves aren't shortened by collapse_short_branch_lengths()
    # so shorten them manually by iterating over all leaves
    for node in t.traverse_leaves():
        if node.get_edge_length() <= args.min_branch_length:
            node.set_edge_length(0)

    # Convert branch lengths to substitution scale
    # instead of substitutions-per-site scale
    t.scale_edges(args.num_sites)

    # If we use ascending=True, this breaks TreeCluster clustering
    # (not sure why?)
    t.order("num_descendants_then_edge_length_then_label", ascending=False)
    t.resolve_polytomies()
    print(t)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tree')
    parser.add_argument('--min_branch_length', type=float, default=0.0000021)
    parser.add_argument('--num_sites', type=int, default=29903)
    args = parser.parse_args()
    main(args)
