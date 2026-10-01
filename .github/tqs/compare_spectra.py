from pathlib import Path

def load(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])
    return rows

gr = load("output/tqs_lcdm_cl.dat")
tq = load("output/tqs_model_cl.dat")
n = min(len(gr), len(tq))

with open("output/tqs_spectrum_comparison.csv", "w") as f:
    f.write("ell,TT_pct,EE_pct,TE_pct,phiphi_pct\n")
    for i in range(n):
        def pct(j):
            if j >= len(gr[i]) or j >= len(tq[i]) or gr[i][j] == 0.0:
                return float("nan")
            return 100.0 * (tq[i][j] / gr[i][j] - 1.0)

        f.write(
            f"{int(gr[i][0])},{pct(1)},{pct(2)},{pct(3)},{pct(5)}\n"
        )

print(f"Compared {n} multipoles")
