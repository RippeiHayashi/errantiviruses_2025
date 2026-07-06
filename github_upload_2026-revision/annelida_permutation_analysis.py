#!/usr/bin/env python3

"""
Permutation analysis of host-taxonomic concordance in an Annelida POL RT clade.

The script keeps the POL RT tree topology fixed, identifies highly supported
small internal clades, and randomly shuffles host-species labels across tips.
Because each species maps to a family, this preserves the observed number of
tips per species and per family while testing whether taxonomic clustering is
stronger than expected by chance.

No third-party Python packages are required.
"""

import argparse
import csv
import random
import statistics
from dataclasses import dataclass, field
from typing import Dict, List, Optional


# -----------------------------------------------------------------------------
# Host taxonomy mapping
# -----------------------------------------------------------------------------
#
# Expected tree-tip format:
#   Annelida_Polychaeta_Pdum_errantivirus_1
#   Annelida_Sipuncula_Snud_errantivirus_1
#
# The third underscore-delimited field is the host code. Add entries here if
# additional Annelida hosts are included in a future version of the tree.

HOST_CODE_TO_TAXON = {
    "Apac": {
        "species": "Amphiduros_pacificus",
        "family": "Hesionidae",
    },
    "Avir": {
        "species": "Alitta_virens",
        "family": "Nereididae",
    },
    "Himp": {
        "species": "Harmothoe_impar",
        "family": "Polynoidae",
    },
    "Lcla": {
        "species": "Lepidonotus_clava",
        "family": "Polynoidae",
    },
    "Pdum": {
        "species": "Platynereis_dumerilii",
        "family": "Nereididae",
    },
    "Slim": {
        "species": "Sthenelais_limicola",
        "family": "Sigalionidae",
    },
    "Snud": {
        "species": "Sipunculus_nudus",
        "family": "Sipunculidae",
    },
}


# -----------------------------------------------------------------------------
# Minimal Newick parser
# Handles IQ-TREE-style comments such as:
#   ):0.0123[100]
# -----------------------------------------------------------------------------

@dataclass
class Node:
    name: str = ""
    length: Optional[float] = None
    support: Optional[float] = None
    children: List["Node"] = field(default_factory=list)
    tips: List[str] = field(default_factory=list)

    def is_leaf(self) -> bool:
        return len(self.children) == 0


class NewickParser:
    def __init__(self, text: str):
        self.s = text.strip()
        self.i = 0

    def peek(self) -> str:
        return self.s[self.i] if self.i < len(self.s) else ""

    def consume_ws(self) -> None:
        while self.i < len(self.s) and self.s[self.i].isspace():
            self.i += 1

    def parse(self) -> Node:
        node = self.parse_subtree()
        self.consume_ws()
        if self.peek() == ";":
            self.i += 1
        self.consume_ws()
        if self.i != len(self.s):
            raise ValueError(f"Unexpected trailing text at position {self.i}")
        return node

    def parse_subtree(self) -> Node:
        self.consume_ws()

        if self.peek() == "(":
            self.i += 1
            children: List[Node] = []

            while True:
                children.append(self.parse_subtree())
                self.consume_ws()

                if self.peek() == ",":
                    self.i += 1
                    continue
                if self.peek() == ")":
                    self.i += 1
                    break
                raise ValueError(f"Expected ',' or ')' at position {self.i}")

            node = Node(children=children)
            self.parse_suffix(node, internal=True)
            return node

        node = Node()
        start = self.i
        while self.i < len(self.s) and self.s[self.i] not in ":,);[":
            self.i += 1
        node.name = self.s[start:self.i].strip()
        if not node.name:
            raise ValueError(f"Missing leaf label at position {start}")
        self.parse_suffix(node, internal=False)
        return node

    def parse_suffix(self, node: Node, internal: bool = False) -> None:
        self.consume_ws()

        # Optional internal node label before the branch length.
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

        # Optional branch length.
        if self.peek() == ":":
            self.i += 1
            start = self.i
            while self.i < len(self.s) and self.s[self.i] not in "[,);":
                self.i += 1
            value = self.s[start:self.i].strip()
            if value:
                try:
                    node.length = float(value)
                except ValueError:
                    pass

        self.consume_ws()

        # Optional IQ-TREE/Newick comments, e.g. [100].
        # Numeric comments on internal nodes are interpreted as bootstrap support.
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


# -----------------------------------------------------------------------------
# Tree utilities
# -----------------------------------------------------------------------------

def annotate_tips(node: Node) -> List[str]:
    if node.is_leaf():
        node.tips = [node.name]
    else:
        tips: List[str] = []
        for child in node.children:
            annotate_tips(child)
            tips.extend(child.tips)
        node.tips = tips
    return node.tips


def iter_internal_nodes(node: Node):
    if not node.is_leaf():
        yield node
        for child in node.children:
            yield from iter_internal_nodes(child)


def host_code_from_tip(tip_label: str) -> str:
    fields = tip_label.split("_")
    if len(fields) < 4 or fields[0] != "Annelida":
        raise ValueError(f"Cannot infer Annelida host code from tip label: {tip_label}")
    return fields[2]


def count_items(values: List[str]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return counts


def format_counts(counter: Dict[str, int]) -> str:
    return "; ".join(f"{key}:{value}" for key, value in sorted(counter.items()))


def quantile(values: List[int], q: float):
    values = sorted(values)
    if not values:
        return None
    pos = (len(values) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(values) - 1)
    frac = pos - lo
    return values[lo] * (1 - frac) + values[hi] * frac


# -----------------------------------------------------------------------------
# Scoring
# -----------------------------------------------------------------------------

def score_clades(clades: List[Node], host_code_by_tip: Dict[str, str]):
    totals = {
        "same_species": 0,
        "same_family": 0,
        "dominant_species_80pct": 0,
        "dominant_family_80pct": 0,
    }

    rows = []

    for node_id, node in enumerate(clades, start=1):
        tips = node.tips
        host_codes = [host_code_by_tip[tip] for tip in tips]
        species = [HOST_CODE_TO_TAXON[code]["species"] for code in host_codes]
        families = [HOST_CODE_TO_TAXON[code]["family"] for code in host_codes]

        host_code_counts = count_items(host_codes)
        species_counts = count_items(species)
        family_counts = count_items(families)

        n_tips = len(tips)
        max_species_fraction = max(species_counts.values()) / n_tips
        max_family_fraction = max(family_counts.values()) / n_tips

        same_species = len(species_counts) == 1
        same_family = len(family_counts) == 1
        dominant_species_80pct = max_species_fraction >= 0.80
        dominant_family_80pct = max_family_fraction >= 0.80

        totals["same_species"] += int(same_species)
        totals["same_family"] += int(same_family)
        totals["dominant_species_80pct"] += int(dominant_species_80pct)
        totals["dominant_family_80pct"] += int(dominant_family_80pct)

        rows.append({
            "node_id": node_id,
            "bootstrap": node.support,
            "n_tips": n_tips,
            "host_code_counts": format_counts(host_code_counts),
            "species_counts": format_counts(species_counts),
            "family_counts": format_counts(family_counts),
            "same_species": same_species,
            "same_family": same_family,
            "dominant_species_80pct": dominant_species_80pct,
            "dominant_species_fraction": round(max_species_fraction, 4),
            "dominant_family_80pct": dominant_family_80pct,
            "dominant_family_fraction": round(max_family_fraction, 4),
            "tip_labels": ";".join(tips),
        })

    return totals, rows


# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Permutation analysis of host-taxonomic concordance in an "
            "Annelida errantivirus POL RT clade."
        )
    )
    parser.add_argument("treefile", help="Newick tree file")
    parser.add_argument("--bootstrap", type=float, default=95.0)
    parser.add_argument("--min-tips", type=int, default=3)
    parser.add_argument("--max-tips", type=int, default=10)
    parser.add_argument("--nperm", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument(
        "--observed-out",
        default="observed_supported_clades_annelida.tsv",
    )
    parser.add_argument(
        "--summary-out",
        default="permutation_summary_annelida.tsv",
    )
    args = parser.parse_args()

    with open(args.treefile, encoding="utf-8") as handle:
        newick = handle.read().strip()

    root = NewickParser(newick).parse()
    annotate_tips(root)

    all_tips = root.tips
    host_code_by_tip = {tip: host_code_from_tip(tip) for tip in all_tips}

    missing_codes = sorted(set(host_code_by_tip.values()) - set(HOST_CODE_TO_TAXON))
    if missing_codes:
        raise ValueError(
            "These host codes are missing from HOST_CODE_TO_TAXON: "
            + ", ".join(missing_codes)
            + ". Add them to the mapping near the top of the script."
        )

    # Select supported small internal clades.
    clades: List[Node] = []
    for node in iter_internal_nodes(root):
        if len(node.children) < 2:
            continue
        if node.support is None:
            continue
        if node.support < args.bootstrap:
            continue
        if args.min_tips <= len(node.tips) <= args.max_tips:
            clades.append(node)

    if not clades:
        raise ValueError(
            "No internal clades matched the selected bootstrap and tip-count thresholds."
        )

    observed_counts, observed_rows = score_clades(clades, host_code_by_tip)

    # Write observed clades.
    with open(args.observed_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(observed_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(observed_rows)

    # Permutation test: shuffle species-code labels across the fixed tree.
    # This preserves the exact observed number of tips per host species.
    rng = random.Random(args.seed)
    tip_order = list(all_tips)
    observed_host_codes = [host_code_by_tip[tip] for tip in tip_order]

    null = {
        "same_species": [],
        "same_family": [],
        "dominant_species_80pct": [],
        "dominant_family_80pct": [],
    }

    for _ in range(args.nperm):
        shuffled_codes = observed_host_codes[:]
        rng.shuffle(shuffled_codes)
        permuted_host_code_by_tip = dict(zip(tip_order, shuffled_codes))

        permuted_counts, _ = score_clades(clades, permuted_host_code_by_tip)
        for metric in null:
            null[metric].append(permuted_counts[metric])

    # Summarize permutation results.
    summary_rows = []
    for metric, observed in observed_counts.items():
        values = null[metric]
        p_empirical = (1 + sum(value >= observed for value in values)) / (
            args.nperm + 1
        )

        summary_rows.append({
            "metric": metric,
            "observed": observed,
            "n_scored_clades": len(clades),
            "null_mean": round(statistics.mean(values), 5),
            "null_median": quantile(values, 0.50),
            "null_2.5pct": quantile(values, 0.025),
            "null_97.5pct": quantile(values, 0.975),
            "empirical_p_ge_observed": p_empirical,
            "n_permutations": args.nperm,
            "random_seed": args.seed,
        })

    with open(args.summary_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    # Console report.
    species = {HOST_CODE_TO_TAXON[code]["species"] for code in host_code_by_tip.values()}
    families = {HOST_CODE_TO_TAXON[code]["family"] for code in host_code_by_tip.values()}

    print(f"Tips: {len(all_tips)}")
    print(f"Host species: {len(species)}")
    print(f"Host families: {len(families)}")
    print(f"Scored internal clades: {len(clades)}")
    print()
    print("Observed:")
    for metric, value in observed_counts.items():
        print(f"  {metric}: {value}")
    print()
    print(f"Wrote: {args.observed_out}")
    print(f"Wrote: {args.summary_out}")


if __name__ == "__main__":
    main()
