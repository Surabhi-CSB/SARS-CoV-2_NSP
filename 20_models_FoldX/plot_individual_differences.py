#!/usr/bin/env python3
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ─── USER SETTINGS ───────────────────────────────────────────────────────────
ROOT = Path(".")
STRUCTURES = ["7jlt", "2ahm", "7bv1", "6yyt"]

# define map for 7jlt
jlt_map = {
    1: "CB114A",   2: "CD114A",   3: "CB114A,CD114A",
    4: "CB114Y",   5: "CD114Y",   6: "CB114Y,CD114Y",
    7: "CB114R",   8: "CD114R",   9: "CB114R,CD114R",
   10: "CB114S",  11: "CD114S",  12: "CB114S,CD114S",
}

# define map for 2ahm
ahm_map = {
    i: name for i, name in enumerate([
        "CE114A","CF114A","CG114A","CH114A","CEH114A",
        "CE114Y","CF114Y","CG114Y","CH114Y","CEH114Y",
        "CE114R","CF114R","CG114R","CH114R","CEH114R",
        "CE114S","CF114S","CG114S","CH114S","CEH114S",
    ], start=1)
}

# now assemble the full MUT_MAP
MUT_MAP = {
    "7jlt": jlt_map,
    "2ahm": ahm_map,
    # 7bv1 and 6yyt use the same chain‐B/D mapping as 7jlt
    "7bv1": jlt_map,
    "6yyt": jlt_map,
}

# which energy terms to plot
ENERGY_TERMS = [
    "total energy", "Van der Waals", "Electrostatics",
    "Solvation Polar", "Solvation Hydrophobic"
]

# ─── FUNCTION TO READ A FXOUT ────────────────────────────────────────────────
def read_fxout(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", skiprows=8, engine="python")
    # drop empty columns if any
    df = df.loc[:, df.columns.notna()]
    return df

# ─── MAIN LOOP: parse + plot per structure ─────────────────────────────────
for struct in STRUCTURES:
    fx = ROOT/struct/f"Dif_{struct}_Repair.fxout"
    if not fx.exists():
        print(f"⚠️  {fx} not found, skipping {struct}")
        continue

    df = read_fxout(fx)
    # extract repair index from Pdb name
    pat = re.compile(fr"{struct}_Repair_(\d+)_\d+\.pdb")
    df["RepairIdx"] = df["Pdb"].str.extract(pat).astype(int)
    df["Mutation"]  = df["RepairIdx"].map(MUT_MAP[struct])

    # melt to long form
    long = (
        df[["Pdb","Mutation","RepairIdx"] + ENERGY_TERMS]
        .melt(
            id_vars=["Pdb","Mutation","RepairIdx"],
            value_vars=ENERGY_TERMS,
            var_name="Term",
            value_name="ΔΔG"
        )
    )

    # plotting
    plt.figure(figsize=(12,6))
    ax = plt.gca()
    terms = ENERGY_TERMS
    x = np.arange(len(terms))
    width = 0.15

    muts = sorted(long["Mutation"].unique())
    colors = plt.cm.tab10.colors

    for i, mut in enumerate(muts):
        sub = long[long["Mutation"]==mut]
        grouped = sub.groupby("Term")["ΔΔG"].apply(list)
        for ti, term in enumerate(terms):
            vals = grouped.get(term, [])
            jitter = (np.random.rand(len(vals)) - 0.5)*width
            ax.scatter(
                x[ti] + (i - len(muts)/2)*width + jitter,
                vals,
                label=mut if ti==0 else "",
                color=colors[i % len(colors)],
                alpha=0.6, edgecolor="k", s=30
            )

    ax.set_xticks(x)
    ax.set_xticklabels(terms, rotation=45, ha="right")
    ax.set_ylabel("ΔΔG (kcal·mol⁻¹)")
    ax.set_title(f"{struct.upper()}: individual replicate ΔΔG")
    ax.legend(ncol=2, bbox_to_anchor=(1.02,1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(f"{struct}_individual_ΔΔG.png", dpi=300)
    plt.show()

