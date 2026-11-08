#!/usr/bin/env python3
"""Task 4 - Statistics and visualization of HBV mutations.

Reads the Task 3 mutation table and writes summary CSV/JSON files and figures.

Example:
    python3 statistics.py
    python3 statistics.py --mutations results/mutations.csv --top 30
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no display needed (servers, cron, CI)
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

log = logging.getLogger("statistics")

# Internal name -> accepted column names in the input CSV.
ALIASES = {
    "accession": ["sequence_accession", "accession", "sample_id"],
    "position": ["position", "reference_position"],
    "ref": ["reference_nucleotide", "ref", "reference_base"],
    "alt": ["sample_nucleotide", "observed_nucleotide", "alt", "observed_base"],
    "type": ["mutation_type", "type", "variation_type"],
}
AMBIGUOUS = set("RYKMSWBDHVN")
GENOTYPE_COLUMNS = ["genotype", "Genotype", "hbv_genotype"]


def standardize(df):
    """Rename input columns to internal names; fail clearly if one is missing."""
    rename, missing = {}, []
    for internal, names in ALIASES.items():
        found = next((n for n in names if n in df.columns), None)
        if found is None:
            missing.append(f"{internal} (one of {names})")
        else:
            rename[found] = internal
    if missing:
        raise ValueError("Missing columns in mutation table: " + "; ".join(missing))
    out = df.rename(columns=rename)[list(ALIASES)].copy()
    out["position"] = pd.to_numeric(out["position"], errors="coerce")
    out = out.dropna(subset=["position"])
    out["position"] = out["position"].astype(int)
    for col in ("ref", "alt", "type"):
        out[col] = out[col].astype(str).str.strip()
    return out


def split_uncertain(df):
    """Separate confirmed mutations from uncertain (ambiguous/IUPAC) calls."""
    is_unc = (
        df["type"].str.lower().str.contains("uncertain|ambiguous")
        | df["alt"].str.upper().apply(lambda s: bool(set(s) & AMBIGUOUS))
    )
    return df[~is_unc].copy(), df[is_unc].copy()


def add_label(df):
    df = df.copy()
    df["mutation"] = (
        df["position"].astype(str) + ":" + df["ref"] + ">" + df["alt"]
    )
    return df


def compute_frequencies(confirmed, denominator):
    """Frequency = distinct sequences carrying the mutation / denominator."""
    grouped = (
        add_label(confirmed)
        .groupby(["mutation", "position", "ref", "alt", "type"])["accession"]
        .nunique()
        .reset_index(name="n_sequences")
    )
    grouped["frequency"] = grouped["n_sequences"] / denominator
    return grouped.sort_values(
        ["n_sequences", "position"], ascending=[False, True]
    ).reset_index(drop=True)


def type_counts(confirmed):
    return (
        confirmed.groupby("type")
        .agg(n_mutations=("accession", "size"), n_sequences=("accession", "nunique"))
        .reset_index()
        .sort_values("n_mutations", ascending=False)
    )


def find_accession_column(df):
    return next((c for c in ALIASES["accession"] if c in df.columns), None)


def pick_denominator(mut_acc, metadata, skipped, forced):
    """Return (denominator, explanation). Priority: --total-sequences, then
    metadata minus skipped, then distinct accessions in the mutation table."""
    if forced:
        return forced, "user-supplied --total-sequences"
    if metadata is not None:
        col = find_accession_column(metadata)
        if col:
            accs = set(metadata[col].astype(str))
            if skipped is not None:
                scol = find_accession_column(skipped)
                if scol:
                    accs -= set(skipped[scol].astype(str))
            return len(accs), "metadata sequences minus skipped sequences"
    return len(mut_acc), (
        "WARNING: only sequences with >=1 mutation (frequencies overestimated)"
    )


def genotype_table(confirmed, metadata, skipped, top_n):
    """Mutation frequency per genotype. Returns None if genotype is unusable."""
    if metadata is None:
        return None
    acc_col = find_accession_column(metadata)
    geno_col = next((c for c in GENOTYPE_COLUMNS if c in metadata.columns), None)
    if acc_col is None or geno_col is None:
        log.warning("No accession/genotype column in metadata: genotype comparison skipped")
        return None
    meta = metadata[[acc_col, geno_col]].rename(
        columns={acc_col: "accession", geno_col: "genotype"}
    )
    meta["accession"] = meta["accession"].astype(str)
    meta["genotype"] = meta["genotype"].fillna("").astype(str).str.strip()
    meta = meta[meta["genotype"] != ""]
    if skipped is not None and find_accession_column(skipped):
        meta = meta[~meta["accession"].isin(skipped[find_accession_column(skipped)].astype(str))]
    if meta.empty:
        log.warning("Genotype column is empty: genotype comparison skipped")
        return None
    sizes = meta.groupby("genotype")["accession"].nunique()
    merged = add_label(confirmed).assign(accession=lambda d: d["accession"].astype(str))
    merged = merged.merge(meta, on="accession", how="inner")
    top = compute_frequencies(confirmed, 1)["mutation"].head(top_n)
    counts = (
        merged[merged["mutation"].isin(top)]
        .groupby(["mutation", "genotype"])["accession"].nunique()
        .reset_index(name="n_sequences")
    )
    counts["genotype_size"] = counts["genotype"].map(sizes)
    counts["frequency"] = counts["n_sequences"] / counts["genotype_size"]
    return counts


def save_fig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    log.info("Figure written: %s", path)


def plot_top(freq, top_n, path):
    top = freq.head(top_n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.35 * len(top) + 1.5))
    ax.barh(top["mutation"], top["frequency"] * 100)
    ax.set_xlabel("Frequency (% of analysed sequences)")
    ax.set_title(f"Top {len(top)} most frequent mutations")
    save_fig(fig, path)


def plot_genome(confirmed, path, bin_size=50):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.hist(confirmed["position"], bins=range(1, 3300, bin_size))
    ax.set_xlabel("Position on reference NC_003977.2")
    ax.set_ylabel("Number of mutation calls")
    ax.set_title(f"Mutations along the genome (bins of {bin_size} nt)")
    save_fig(fig, path)


def plot_types(types, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(types["type"], types["n_mutations"])
    ax.set_ylabel("Number of mutation calls")
    ax.set_title("Mutations by type")
    save_fig(fig, path)


def plot_heatmap(geno, path):
    pivot = geno.pivot(index="mutation", columns="genotype", values="frequency").fillna(0)
    fig, ax = plt.subplots(figsize=(1 + 0.8 * pivot.shape[1], 0.35 * pivot.shape[0] + 2))
    im = ax.imshow(pivot.values, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]), pivot.columns)
    ax.set_yticks(range(pivot.shape[0]), pivot.index)
    fig.colorbar(im, ax=ax, label="Frequency within genotype")
    ax.set_title("Top mutations by genotype")
    save_fig(fig, path)


def read_csv_or_exit(path, required=True):
    path = Path(path)
    if not path.is_file():
        if required:
            log.error("Input file not found: %s", path)
            sys.exit(2)
        return None
    return pd.read_csv(path, dtype=str)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Task 4: HBV mutation statistics and figures.")
    p.add_argument("--mutations", default="results/mutations.csv")
    p.add_argument("--metadata", default="data/processed/hbv_metadata.csv")
    p.add_argument("--skipped", default="results/mutation_analysis_skipped.csv")
    p.add_argument("--outdir", default="results")
    p.add_argument("--top", type=int, default=20, help="number of top mutations to plot")
    p.add_argument("--total-sequences", type=int, default=None,
                   help="number of analysed sequences (overrides automatic estimate)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    try:
        raw = read_csv_or_exit(args.mutations)
        data = standardize(raw)
    except ValueError as err:
        log.error("%s", err)
        return 2
    metadata = read_csv_or_exit(args.metadata, required=False)
    skipped = read_csv_or_exit(args.skipped, required=False)
    if metadata is None:
        log.warning("Metadata not found (%s): no genotype comparison", args.metadata)

    confirmed, uncertain = split_uncertain(data)
    denom, basis = pick_denominator(
        data["accession"].unique(), metadata, skipped, args.total_sequences
    )
    if denom <= 0:
        log.error("No sequences to compute frequencies on")
        return 2
    log.info("Denominator = %d (%s)", denom, basis)

    freq = compute_frequencies(confirmed, denom)
    types = type_counts(confirmed)
    freq.to_csv(outdir / "summary.csv", index=False)
    types.to_csv(outdir / "mutation_types.csv", index=False)

    if not freq.empty:
        plot_top(freq, args.top, figdir / "top_mutations.png")
        plot_genome(confirmed, figdir / "mutations_along_genome.png")
        plot_types(types, figdir / "mutation_types.png")

    geno = genotype_table(confirmed, metadata, skipped, args.top)
    if geno is not None and not geno.empty:
        geno.to_csv(outdir / "genotype_comparison.csv", index=False)
        plot_heatmap(geno, figdir / "genotype_heatmap.png")

    summary = {
        "denominator": int(denom),
        "denominator_basis": basis,
        "n_confirmed_calls": int(len(confirmed)),
        "n_uncertain_calls_excluded": int(len(uncertain)),
        "n_distinct_mutations": int(len(freq)),
        "genotype_comparison": geno is not None and not geno.empty,
    }
    (outdir / "statistics_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("Done: %s", summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
