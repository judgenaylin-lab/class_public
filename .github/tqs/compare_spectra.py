from pathlib import Path

def load(path):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])
    return rows

def find_cl(prefix):
    out = Path("output")
    exact = out / f"{prefix}cl.dat"
    if exact.exists():
        return exact
    candidates = sorted(out.glob(f"{prefix}*cl*.dat"))
    if not candidates:
        available = ", ".join(p.name for p in sorted(out.glob("*")))
        raise FileNotFoundError(
            f"No CLASS C_l file found for prefix {prefix}. Available: {available}"
        )
    unlensed = [p for p in candidates if "lensed" not in p.name]
    return unlensed[0] if unlensed else candidates[0]

gr_path = find_cl("tqs_lcdm_")
tq_path = find_cl("tqs_model_")
print(f"Using baseline spectra: {gr_path}")
print(f"Using TQS spectra: {tq_path}")

gr = load(gr_path)
tq = load(tq_path)
n = min(len(gr), len(tq))

with open("output/tqs_spectrum_comparison.csv", "w") as f:
    f.write("ell,TT_pct,EE_pct,TE_pct,phiphi_pct\n")
    for i in range(n):
        def pct(j):
            if j >= len(gr[i]) or j >= len(tq[i]) or gr[i][j] == 0.0:
                return float("nan")
            return 100.0 * (tq[i][j] / gr[i][j] - 1.0)
        f.write(f"{int(gr[i][0])},{pct(1)},{pct(2)},{pct(3)},{pct(5)}\n")

print(f"Compared {n} multipoles")
