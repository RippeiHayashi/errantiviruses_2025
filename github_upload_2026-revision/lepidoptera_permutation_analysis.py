#!/usr/bin/env python3

import argparse
import csv
import random
import statistics
from dataclasses import dataclass, field
from typing import List, Optional


# -----------------------------
# Family-to-superfamily mapping
# -----------------------------

FAMILY_TO_SUPERFAMILY = {
    "Adelidae": "Adeloidea",
    "Incurvariidae": "Incurvarioidea",
    "Tineidae": "Tineoidea",

    "Yponomeutidae": "Yponomeutoidea",
    "Ypsolophidae": "Yponomeutoidea",
    "Plutellidae": "Yponomeutoidea",

    "Gelechiidae": "Gelechioidea",
    "Coleophoridae": "Gelechioidea",
    "Blastobasidae": "Gelechioidea",
    "Depressariidae": "Gelechioidea",
    "Oecophoridae": "Gelechioidea",

    "Carposinidae": "Carposinoidea",
    "Tortricidae": "Tortricoidea",
    "Zygaenidae": "Zygaenoidea",
    "Cossidae": "Cossoidea",
    "Sesiidae": "Sesioidea",

    "Pterophoridae": "Pterophoroidea",

    "Pyralidae": "Pyraloidea",
    "Crambidae": "Pyraloidea",

    "Papilionidae": "Papilionoidea",
    "Pieridae": "Papilionoidea",
    "Lycaenidae": "Papilionoidea",
    "Nymphalidae": "Papilionoidea",
    "Hesperiidae": "Papilionoidea",

    "Drepanidae": "Drepanoidea",
    "Geometridae": "Geometroidea",
    "Lasiocampidae": "Lasiocampoidea",

    "Bombycidae": "Bombycoidea",
    "Saturniidae": "Bombycoidea",
    "Sphingidae": "Bombycoidea",

    "Notodontidae": "Noctuoidea",
    "Erebidae": "Noctuoidea",
    "Nolidae": "Noctuoidea",
    "Noctuidae": "Noctuoidea",
}


# -----------------------------
# Minimal Newick parser
# Handles IQ-TREE-style comments:
# ):0.0123[100]
# -----------------------------

@dataclass
class Node:
    name: str = ""
    length: Optional[float] = None
    support: Optional[float] = None
    children: List["Node"] = field(default_factory=list)
    tips: List[str] = field(default_factory=list)

    def is_leaf(self):
        return len(self.children) == 0


class NewickParser:
    def __init__(self, text):
        self.s = text.strip()
        self.i = 0

    def peek(self):
        return self.s[self.i] if self.i < len(self.s) else ""

    def consume_ws(self):
        while self.i < len(self.s) and self.s[self.i].isspace():
            self.i += 1

    def parse(self):
        node = self.parse_subtree()
        self.consume_ws()
        if self.peek() == ";":
            self.i += 1
        return node

    def parse_subtree(self):
        self.consume_ws()

        if self.peek() == "(":
            self.i += 1
            children = []

            while True:
                children.append(self.parse_subtree())
                self.consume_ws()

                if self.peek() == ",":
                    self.i += 1
                    continue
                elif self.peek() == ")":
                    self.i += 1
                    break
                else:
                    raise ValueError(f"Expected ',' or ')' at position {self.i}")

            node = Node(children=children)
            self.parse_suffix(node, internal=True)
            return node

        else:
            node = Node()
            start = self.i
            while self.i < len(self.s) and self.s[self.i] not in ":,);[":
                self.i += 1
            node.name = self.s[start:self.i].strip()
            self.parse_suffix(node, internal=False)
            return node

    def parse_suffix(self, node, internal=False):
        self.consume_ws()

        # Optional internal node label before branch length.
        if internal and self.peek() not in ":,);[":
            start = self.i
            while self.i < len(self.s) and self.s[self.i] not in ":,);[":
                self.i += 1
            label = self.s[start:self.i].strip()
            if label:
                node.name = label
                try:
                    node.support = float(label)
                except ValueError:
                    pass

        self.consume_ws()

        # Branch length.
        if self.peek() == ":":
            self.i += 1
            start = self.i
            while self.i < len(self.s) and self.s[self.i] not in "[,);":
                self.i += 1
            val = self.s[start:self.i].strip()
            if val:
                try:
                    node.length = float(val)
                except ValueError:
                    pass

        self.consume_ws()

        # IQ-TREE/Newick comments, e.g. [100].
        # For internal nodes, numeric comments are treated as bootstrap support.
        while self.peek() == "[":
            self.i += 1
            start = self.i
            while self.i < len(self.s) and self.s[self.i] != "]":
                self.i += 1
            comment = self.s[start:self.i].strip()
            if self.peek() == "]":
                self.i += 1

            if internal:
                try:
                    node.support = float(comment)
                except ValueError:
                    pass

            self.consume_ws()


# -----------------------------
# Tree utilities
# -----------------------------

def annotate_tips(node):
    if node.is_leaf():
        node.tips = [node.name]
    else:
        tips = []
        for child in node.children:
            annotate_tips(child)
            tips.extend(child.tips)
        node.tips = tips
    return node.tips


def iter_internal_nodes(node):
    if not node.is_leaf():
        yield node
        for child in node.children:
            yield from iter_internal_nodes(child)


def family_from_tip(tip_label):
    # Expected format:
    # Lepidoptera_Family_GenusSpecies_errantivirus_N
    fields = tip_label.split("_")
    if len(fields) < 2 or fields[0] != "Lepidoptera":
        raise ValueError(f"Cannot infer family from tip label: {tip_label}")
    return fields[1]


def count_items(values):
    out = {}
    for v in values:
        out[v] = out.get(v, 0) + 1
    return out


def format_counts(counter):
    return "; ".join(f"{k}:{v}" for k, v in sorted(counter.items()))


# -----------------------------
# Scoring
# -----------------------------

def score_clades(clades, family_by_tip):
    totals = {
        "same_family": 0,
        "same_superfamily": 0,
        "dominant_superfamily_80pct": 0,
    }

    rows = []

    for i, node in enumerate(clades, start=1):
        tips = node.tips
        families = [family_by_tip[t] for t in tips]
        superfamilies = [FAMILY_TO_SUPERFAMILY[f] for f in families]

        family_counts = count_items(families)
        superfamily_counts = count_items(superfamilies)

        n = len(tips)
        max_sf_count = max(superfamily_counts.values())
        max_sf_fraction = max_sf_count / n

        same_family = len(family_counts) == 1
        same_superfamily = len(superfamily_counts) == 1
        dominant_superfamily_80pct = max_sf_fraction >= 0.80

        totals["same_family"] += int(same_family)
        totals["same_superfamily"] += int(same_superfamily)
        totals["dominant_superfamily_80pct"] += int(dominant_superfamily_80pct)

        rows.append({
            "node_id": i,
            "bootstrap": node.support,
            "n_tips": n,
            "family_counts": format_counts(family_counts),
            "superfamily_counts": format_counts(superfamily_counts),
            "same_family": same_family,
            "same_superfamily": same_superfamily,
            "dominant_superfamily_80pct": dominant_superfamily_80pct,
            "dominant_superfamily_fraction": round(max_sf_fraction, 4),
            "tip_labels": ";".join(tips),
        })

    return totals, rows


def quantile(values, q):
    values = sorted(values)
    if not values:
        return None
    pos = (len(values) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(values) - 1)
    frac = pos - lo
    return values[lo] * (1 - frac) + values[hi] * frac


# -----------------------------
# Main analysis
# -----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Permutation analysis of host-taxonomic concordance in Lepidoptera POL RT clade."
    )
    parser.add_argument("treefile", help="Newick tree file")
    parser.add_argument("--bootstrap", type=float, default=95.0)
    parser.add_argument("--min-tips", type=int, default=3)
    parser.add_argument("--max-tips", type=int, default=10)
    parser.add_argument("--nperm", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--observed-out", default="observed_supported_clades.tsv")
    parser.add_argument("--summary-out", default="permutation_summary.tsv")
    args = parser.parse_args()

    with open(args.treefile) as f:
        newick = f.read().strip()

    root = NewickParser(newick).parse()
    annotate_tips(root)

    all_tips = root.tips
    family_by_tip = {tip: family_from_tip(tip) for tip in all_tips}

    missing = sorted(set(family_by_tip.values()) - set(FAMILY_TO_SUPERFAMILY))
    if missing:
        raise ValueError(
            "These families are missing from FAMILY_TO_SUPERFAMILY mapping: "
            + ", ".join(missing)
        )

    # Select supported small internal clades.
    clades = []
    for node in iter_internal_nodes(root):
        if len(node.children) < 2:
            continue
        if node.support is None:
            continue
        if node.support < args.bootstrap:
            continue
        if args.min_tips <= len(node.tips) <= args.max_tips:
            clades.append(node)

    observed_counts, observed_rows = score_clades(clades, family_by_tip)

    # Write observed clades.
    with open(args.observed_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(observed_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(observed_rows)

    # Permutation test.
    rng = random.Random(args.seed)
    tip_order = list(all_tips)
    observed_family_labels = [family_by_tip[t] for t in tip_order]

    null = {
        "same_family": [],
        "same_superfamily": [],
        "dominant_superfamily_80pct": [],
    }

    for _ in range(args.nperm):
        shuffled = observed_family_labels[:]
        rng.shuffle(shuffled)
        permuted_family_by_tip = dict(zip(tip_order, shuffled))

        perm_counts, _ = score_clades(clades, permuted_family_by_tip)

        for metric in null:
            null[metric].append(perm_counts[metric])

    # Summarize permutation results.
    summary_rows = []
    for metric, observed in observed_counts.items():
        vals = null[metric]
        p_empirical = (1 + sum(v >= observed for v in vals)) / (args.nperm + 1)

        summary_rows.append({
            "metric": metric,
            "observed": observed,
            "n_scored_clades": len(clades),
            "null_mean": round(statistics.mean(vals), 5),
            "null_median": quantile(vals, 0.50),
            "null_2.5pct": quantile(vals, 0.025),
            "null_97.5pct": quantile(vals, 0.975),
            "empirical_p_ge_observed": p_empirical,
            "n_permutations": args.nperm,
            "random_seed": args.seed,
        })

    with open(args.summary_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    # Console report.
    print(f"Tips: {len(all_tips)}")
    print(f"Families: {len(set(family_by_tip.values()))}")
    print(f"Scored internal clades: {len(clades)}")
    print()
    print("Observed:")
    for k, v in observed_counts.items():
        print(f"  {k}: {v}")
    print()
    print(f"Wrote: {args.observed_out}")
    print(f"Wrote: {args.summary_out}")


if __name__ == "__main__":
    main()
