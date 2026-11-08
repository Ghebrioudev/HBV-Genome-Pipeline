#!/usr/bin/env python3
"""Task 4 - Statistics and visualization of HBV mutations.

Reads the Task 3 mutation table and writes summary CSV/JSON files and figures.

The table can be very large (the full dataset produces >10 M rows), so the
module streams it in chunks and keeps only aggregate counters in memory.

Example:
    python3 statistics.py
    python3 statistics.py --mutations results/mutations.csv --top 30
"""
import argparse
import json
import logging
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no display needed (servers, cron, CI)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
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
AMBIGUOUS_RE = "[" + "".join(sorted(AMBIGUOUS)) + "]"
GENOTYPE_COLUMNS = ["genotype", "Genotype", "hbv_genotype"]
MAX_REFERENCE_POSITION = 3182
DEFAULT_CHUNKSIZE = 250_000
BIN_SIZE = 50

_ANCHOR_RE = re.compile(r"^\d+\^\d+$")


def _rename_map(df):
    rename, missing = {}, []
    for internal, names in ALIASES.items():
        found = next((n for n in names if n in df.columns), None)
        if found is None:
            missing.append(f"{internal} (one of {names})")
        else:
            rename[found] = internal
    return rename, missing


def standardize(df):
    """Rename input columns to internal names; fail clearly if one is missing.

    Numeric positions become integers. Insertion anchors such as ``125^126``
    are dropped here for backwards compatibility; :func:`standardize_chunk`
    keeps them for the streaming path.
    """
    rename, missing = _rename_map(df)
    if missing:
        raise ValueError("Missing columns in mutation table: " + "; ".join(missing))
    out = df.rename(columns=rename)[list(ALIASES)].copy()
    out["position"] = pd.to_numeric(out["position"], errors="coerce")
    out = out.dropna(subset=["position"])
    out["position"] = out["position"].astype(int)
    for col in ("ref", "alt", "type"):
        out[col] = out[col].astype(str).str.strip()
    return out


def standardize_chunk(df):
    """Like :func:`standardize`, but preserves insertion anchors.

    Adds ``position_key`` (the first reference coordinate) for sorting and
    plotting, and keeps ``position`` as the original label.
    """
    rename, missing = _rename_map(df)
    if missing:
        raise ValueError("Missing columns in mutation table: " + "; ".join(missing))
    out = df.rename(columns=rename)[list(ALIASES)].copy()
    out["position"] = out["position"].astype(str).str.strip()
    keys = pd.to_numeric(out["position"].str.split("^").str[0], errors="coerce")
    keep = keys.notna()
    out = out[keep].copy()
    out["position_key"] = keys[keep].astype("int64")
    for col in ("ref", "alt", "type"):
        out[col] = out[col].astype(str).str.strip()
    return out


def uncertain_mask(df):
    """Vectorized equivalent of the ambiguity test in :func:`split_uncertain`."""
    return df["type"].str.lower().str.contains("uncertain|ambiguous", regex=True) | df[
        "alt"
    ].str.upper().str.contains(AMBIGUOUS_RE, regex=True)


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


def pick_denominator(mut_acc, metadata, skipped, forced, exclude=None):
    """Return (denominator, explanation). Priority: --total-sequences, then
    metadata minus skipped, then distinct accessions in the mutation table."""
    if forced:
        return forced, "user-supplied --total-sequences"
    if metadata is not None:
        col = find_accession_column(metadata)
        if col:
            accs = set(metadata[col].astype(str))
            if exclude:
                accs -= set(exclude)  # the reference is not an analysed sample
            if skipped is not None:
                scol = find_accession_column(skipped)
                if scol:
                    accs -= set(skipped[scol].astype(str))
            return len(accs), "metadata sequences minus skipped sequences"
    return len(mut_acc), (
        "WARNING: only sequences with >=1 mutation (frequencies overestimated)"
    )


def load_metadata_light(path):
    """Load only the accession and (if present) genotype columns."""
    path = Path(path)
    if not path.is_file():
        return None
    header = pd.read_csv(path, nrows=0).columns.tolist()
    acc_col = next((c for c in ALIASES["accession"] if c in header), None)
    if acc_col is None:
        return None
    geno_col = next((c for c in GENOTYPE_COLUMNS if c in header), None)
    columns = [acc_col] + ([geno_col] if geno_col else [])
    return pd.read_csv(path, dtype=str, usecols=columns)


def genotype_table(confirmed, metadata, skipped, top_n):
    """Mutation frequency per genotype. Returns None if genotype is unusable.

    Kept for the in-memory API and tests; :func:`genotype_table_streamed` is the
    chunked equivalent used by :func:`main`.
    """
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
        .groupby(["mutation", "genotype"])["accession"]
        .nunique()
        .reset_index(name="n_sequences")
    )
    counts["genotype_size"] = counts["genotype"].map(sizes)
    counts["frequency"] = counts["n_sequences"] / counts["genotype_size"]
    return counts


def genotype_table_streamed(path, metadata, top_mutations, skipped, chunksize):
    """Chunked version of :func:`genotype_table` for the large table."""
    if metadata is None:
        return None
    acc_col = find_accession_column(metadata)
    geno_col = next((c for c in GENOTYPE_COLUMNS if c in metadata.columns), None)
    if acc_col is None or geno_col is None:
        return None
    meta = metadata[[acc_col, geno_col]].rename(
        columns={acc_col: "accession", geno_col: "genotype"}
    )
    meta["accession"] = meta["accession"].astype(str)
    meta["genotype"] = meta["genotype"].fillna("").astype(str).str.strip()
    meta = meta[meta["genotype"] != ""]
    if meta.empty:
        return None
    sizes = meta.groupby("genotype")["accession"].nunique()
    genotype_of = dict(zip(meta["accession"], meta["genotype"]))
    top = set(top_mutations)
    counts = Counter()
    for chunk in pd.read_csv(path, dtype=str, chunksize=chunksize):
        data = standardize_chunk(chunk)
        conf = data[~uncertain_mask(data)].copy()
        if conf.empty:
            continue
        conf = add_label(conf)
        conf = conf[conf["mutation"].isin(top)]
        if conf.empty:
            continue
        conf["genotype"] = conf["accession"].map(genotype_of)
        conf = conf.dropna(subset=["genotype"])
        conf = conf.drop_duplicates(["accession", "mutation"])
        for (mutation, genotype), n in conf.groupby(["mutation", "genotype"]).size().items():
            counts[(mutation, genotype)] += int(n)
    if not counts:
        return None
    rows = [
        {"mutation": m, "genotype": g, "n_sequences": n}
        for (m, g), n in counts.items()
    ]
    table = pd.DataFrame(rows)
    table["genotype_size"] = table["genotype"].map(sizes)
    table["frequency"] = table["n_sequences"] / table["genotype_size"]
    return table


def stream_aggregate(path, chunksize=DEFAULT_CHUNKSIZE):
    """Scan the mutation table once and return memory-bounded aggregates."""
    edges = np.arange(1, MAX_REFERENCE_POSITION + BIN_SIZE + 1, BIN_SIZE)
    histogram = np.zeros(len(edges) - 1, dtype="int64")
    mutation_counts = {}
    type_rows = Counter()
    type_accessions = defaultdict(set)
    accessions = set()
    reference_accessions = set()
    confirmed_total = 0
    uncertain_total = 0
    for chunk in pd.read_csv(path, dtype=str, chunksize=chunksize):
        if "reference_accession" in chunk.columns:
            reference_accessions.update(chunk["reference_accession"].dropna().astype(str))
        data = standardize_chunk(chunk)
        accessions.update(data["accession"].astype(str))
        mask = uncertain_mask(data)
        uncertain = data[mask]
        confirmed = data[~mask]
        uncertain_total += len(uncertain)
        confirmed_total += len(confirmed)
        if confirmed.empty:
            continue
        confirmed = add_label(confirmed)
        for name, size in confirmed.groupby("type").size().items():
            type_rows[str(name)] += int(size)
        for name, group in confirmed.groupby("type"):
            type_accessions[str(name)].update(group["accession"].astype(str))
        histogram += np.histogram(confirmed["position_key"].to_numpy(), bins=edges)[0]
        # One row per (sequence, mutation) is expected; keep it explicit so a
        # duplicated call cannot inflate the count.
        unique = confirmed.drop_duplicates(["accession", "mutation"])
        grouped = unique.groupby(["position", "ref", "alt", "type"]).size()
        for key, size in grouped.items():
            mutation_counts[key] = mutation_counts.get(key, 0) + int(size)
    return {
        "histogram": histogram,
        "edges": edges,
        "mutation_counts": mutation_counts,
        "type_rows": type_rows,
        "type_accessions": type_accessions,
        "accessions": accessions,
        "reference_accessions": reference_accessions,
        "confirmed_total": confirmed_total,
        "uncertain_total": uncertain_total,
    }


def frequencies_from_counts(mutation_counts, denominator):
    """Turn streaming counters into the Task 4 summary table."""
    columns = ["mutation", "position", "ref", "alt", "type", "n_sequences"]
    if not mutation_counts:
        frame = pd.DataFrame(columns=columns)
        frame["frequency"] = pd.Series(dtype=float)
        return frame
    rows = [
        {
            "mutation": f"{position}:{ref}>{alt}",
            "position": position,
            "ref": ref,
            "alt": alt,
            "type": mutation_type,
            "n_sequences": count,
        }
        for (position, ref, alt, mutation_type), count in mutation_counts.items()
    ]
    frame = pd.DataFrame(rows, columns=columns)
    frame["frequency"] = frame["n_sequences"] / denominator
    sort_key = frame["position"].str.split("^").str[0].astype(int)
    frame = frame.assign(_position_key=sort_key)
    frame = frame.sort_values(
        ["n_sequences", "_position_key"], ascending=[False, True]
    ).drop(columns="_position_key")
    return frame.reset_index(drop=True)


def types_from_counts(type_rows, type_accessions):
    rows = [
        {
            "type": name,
            "n_mutations": int(count),
            "n_sequences": len(type_accessions.get(name, ())),
        }
        for name, count in type_rows.items()
    ]
    frame = pd.DataFrame(rows, columns=["type", "n_mutations", "n_sequences"])
    return frame.sort_values("n_mutations", ascending=False).reset_index(drop=True)


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


def plot_genome_histogram(histogram, edges, path, bin_size=BIN_SIZE):
    centers = (edges[:-1] + edges[1:]) / 2
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(centers, histogram, width=bin_size, align="center")
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
    p.add_argument("--chunksize", type=int, default=DEFAULT_CHUNKSIZE,
                   help="rows per streaming chunk (default 250000)")
    p.add_argument("--total-sequences", type=int, default=None,
                   help="number of analysed sequences (overrides automatic estimate)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    mutations_path = Path(args.mutations)
    if not mutations_path.is_file():
        log.error("Input file not found: %s", mutations_path)
        return 2

    metadata = load_metadata_light(args.metadata)
    skipped = read_csv_or_exit(args.skipped, required=False)
    if metadata is None:
        log.warning("Metadata not found (%s): no genotype comparison", args.metadata)

    try:
        agg = stream_aggregate(mutations_path, args.chunksize)
    except ValueError as err:
        log.error("%s", err)
        return 2

    denom, basis = pick_denominator(
        agg["accessions"], metadata, skipped, args.total_sequences,
        exclude=agg["reference_accessions"],
    )
    if denom <= 0:
        log.error("No sequences to compute frequencies on")
        return 2
    log.info("Denominator = %d (%s)", denom, basis)

    freq = frequencies_from_counts(agg["mutation_counts"], denom)
    types = types_from_counts(agg["type_rows"], agg["type_accessions"])
    freq.to_csv(outdir / "summary.csv", index=False)
    types.to_csv(outdir / "mutation_types.csv", index=False)

    if not freq.empty:
        plot_top(freq, args.top, figdir / "top_mutations.png")
        plot_genome_histogram(agg["histogram"], agg["edges"],
                              figdir / "mutations_along_genome.png")
        plot_types(types, figdir / "mutation_types.png")

    geno = genotype_table_streamed(
        mutations_path, metadata, freq["mutation"].head(args.top),
        skipped, args.chunksize,
    )
    if geno is not None and not geno.empty:
        geno.to_csv(outdir / "genotype_comparison.csv", index=False)
        plot_heatmap(geno, figdir / "genotype_heatmap.png")
    else:
        log.warning("Genotype comparison skipped: metadata has no usable genotype column")

    summary = {
        "denominator": int(denom),
        "denominator_basis": basis,
        "n_confirmed_calls": int(agg["confirmed_total"]),
        "n_uncertain_calls_excluded": int(agg["uncertain_total"]),
        "n_distinct_mutations": int(len(freq)),
        "genotype_comparison": geno is not None and not geno.empty,
    }
    (outdir / "statistics_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("Done: %s", summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
