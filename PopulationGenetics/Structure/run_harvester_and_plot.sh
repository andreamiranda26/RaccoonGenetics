#!/usr/bin/env bash
# run_harvester_and_plot.sh
# One-stop STRUCTURE post-processing:
# - Clone/update structureHarvester
# - Run Harvester (Evanno + CLUMPP files)
# - Report best K
# - Plot meanQ + representative run for each K directly from .out files

source $HOME/miniforge3/bin/activate structviz

set -euo pipefail

### -----------------------------
### CONFIG (change if you prefer)
### -----------------------------

# You can pass --root=/path OR export ROOT_DIR before running.
ROOT_DIR_DEFAULT="/scratch/Raccoon_andrea/gabe_trial/"

# Where your STRUCTURE .out files live:
STRUCT_DIR_REL="analysis/structure"

# Where to put Harvester summaries & plots:
HARVEST_SUBDIR="harvester"
PLOTS_SUBDIR="plots"

# Where to place the harvester code (inside ROOT_DIR by default):
TOOLS_SUBDIR="tools"

### -----------------------------
### ARG PARSING
### -----------------------------
ROOT_DIR="${ROOT_DIR:-$ROOT_DIR_DEFAULT}"
for arg in "$@"; do
  case "$arg" in
    --root=*) ROOT_DIR="${arg#*=}";;
    *) echo "[warn] Unknown arg: $arg";;
  esac
done

STRUCT_DIR="${ROOT_DIR}/${STRUCT_DIR_REL}"
TOOLS_DIR="${ROOT_DIR}/${TOOLS_SUBDIR}"
HARVEST_DIR="${STRUCT_DIR}/${HARVEST_SUBDIR}"
PLOTS_DIR="${STRUCT_DIR}/${PLOTS_SUBDIR}"

mkdir -p "$TOOLS_DIR" "$HARVEST_DIR" "$PLOTS_DIR"

echo "[info] ROOT_DIR        = $ROOT_DIR"
echo "[info] STRUCT_DIR      = $STRUCT_DIR"
echo "[info] HARVEST_DIR     = $HARVEST_DIR"
echo "[info] PLOTS_DIR       = $PLOTS_DIR"

### -----------------------------
### DEP CHECKS
### -----------------------------
need() {
  command -v "$1" >/dev/null 2>&1 || { echo "[error] '$1' not found in PATH"; exit 1; }
}
need git
need python3

# Check python libs for plotting (numpy + matplotlib)
if ! python3 - <<'PY' >/dev/null 2>&1
import numpy, matplotlib.pyplot
PY
then
  echo "[error] Python needs numpy and matplotlib for plotting."
  echo "        Try: pip install --user numpy matplotlib"
  exit 1
fi

### -----------------------------
### GET/UPDATE structureHarvester
### -----------------------------
cd "$TOOLS_DIR"
if [ ! -d "structureHarvester" ]; then
  echo "[info] Cloning structureHarvester..."
  git clone https://github.com/dentearl/structureHarvester.git
else
  echo "[info] Updating structureHarvester..."
  (cd structureHarvester && git pull --ff-only || true)
fi
chmod 755 structureHarvester/structureHarvester.py

### -----------------------------
### RUN HARVESTER
### -----------------------------
echo "[info] Running structureHarvester..."
python3 "structureHarvester/structureHarvester.py" \
  --dir="$STRUCT_DIR" \
  --out="$HARVEST_DIR" \
  --evanno \
  --clumpp

echo "[ok] Harvester finished."
EVANNO_FILE="$HARVEST_DIR/evanno.txt"
if [ -s "$EVANNO_FILE" ]; then
  echo
  echo "========== Evanno (ΔK) summary =========="
  column -t < "$EVANNO_FILE" | sed 's/\t/  /g' | head -n 50
  echo "========================================="
  # Pick best K (max DeltaK over K>1)
  BEST_K=$(awk 'NR>1 && $1>1 { if ($7>max) {max=$7; best=$1} } END {print best+0}' "$EVANNO_FILE")
  if [ "$BEST_K" != "0" ]; then
    echo "[info] Best K by Evanno (ΔK) appears to be: K=${BEST_K}"
  else
    echo "[warn] Could not determine best K (ΔK); check $EVANNO_FILE"
  fi
else
  echo "[warn] Evanno file not found or empty: $EVANNO_FILE"
fi

### -----------------------------
### WRITE A COMPACT PLOTTER
### -----------------------------
PLOTTER="$STRUCT_DIR/plot_structure.py"
cat > "$PLOTTER" <<'PY'
#!/usr/bin/env python3
# Minimal STRUCTURE plotter:
# - finds Raccoon_K{K}_R*.out
# - parses Q from "Estimated membership..." table
# - aligns replicates per K by brute-force permutation
# - plots mean Q and a representative run

import re, os, glob, math, itertools
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT_DIR = os.environ.get("ROOT_DIR", "")
STRUCT_DIR = Path(ROOT_DIR) / "analysis" / "structure" if ROOT_DIR else Path.cwd()
PLOTS_DIR = STRUCT_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

def parse_q(path):
    txt = open(path, 'r', errors='ignore').read()
    m = re.search(r"Estimated membership of individuals in each of K populations", txt)
    if not m: raise ValueError(f"No membership block in {path}")
    lines = txt[m.end():].splitlines()
    rows=[]
    for line in lines:
        line=line.strip()
        if not line:
            if rows: break
            else: continue
        mm = re.match(r"^\d+\s*:\s*(.+)$", line)
        if mm:
            vals=[]
            for tok in mm.group(1).split():
                try: vals.append(float(tok))
                except: break
            if vals: rows.append(vals)
        else:
            if rows: break
    Q = np.array(rows, float)
    Q = Q / Q.sum(axis=1, keepdims=True)
    return Q

def best_perm(Q, ref):
    K = Q.shape[1]
    best=None; bests=-1e9
    for p in itertools.permutations(range(K)):
        s=0.0
        for k in range(K):
            a=Q[:,p[k]]; b=ref[:,k]
            if a.var()==0 or b.var()==0: r=0.0
            else:
                r=np.corrcoef(a,b)[0,1]
                if np.isnan(r): r=0.0
            s+=r
        if s>bests: best, bests = p, s
    return best

def sort_inds(Q):
    maxk=np.argmax(Q, axis=1); maxv=np.max(Q, axis=1)
    order=np.lexsort((-maxv, maxk))
    return Q[order,:], order

def barplot(Q, title, labels, outpng, outpdf):
    N,K = Q.shape
    x=np.arange(N); bottoms=np.zeros(N)
    fig = plt.figure(figsize=(max(8, N/12), 4.0))
    for k in range(K):
        plt.bar(x, Q[:,k], bottom=bottoms, width=0.95, edgecolor='none')
        bottoms+=Q[:,k]
    plt.title(title); plt.ylabel("Ancestry proportion"); plt.ylim(0,1)
    if N<=60:
        plt.xticks(x, labels, rotation=90, fontsize=8)
    else:
        step=math.ceil(N/60); lab=[str(i+1) if i%step==0 else "" for i in range(N)]
        plt.xticks(x, lab, rotation=90, fontsize=7)
    plt.tight_layout(); plt.savefig(outpng, dpi=300); plt.savefig(outpdf); plt.close(fig)

def main():
    # infer Ks from filenames
    files = glob.glob(str(STRUCT_DIR / "Raccoon_K*_R*.out"))
    if not files:
        raise SystemExit(f"[error] No files like Raccoon_K*_R*.out in {STRUCT_DIR}")
    Ks = sorted({int(Path(f).stem.split("_K")[1].split("_")[0]) for f in files})
    for K in Ks:
        fks = sorted(glob.glob(str(STRUCT_DIR / f"Raccoon_K{K}_R*.out")))
        qs=[]
        keep=[]
        for fp in fks:
            try:
                Q=parse_q(fp)
                if Q.shape[1]!=K: continue
                qs.append(Q); keep.append(fp)
            except Exception as e:
                print(f"[warn] {fp}: {e}")
        if not qs: 
            print(f"[warn] K={K}: nothing parsable"); 
            continue
        Ns={q.shape[0] for q in qs}
        if len(Ns)!=1: raise SystemExit(f"[error] Inconsistent N for K={K}: {Ns}")
        # align to the first
        ref=qs[0]; aligned=[ref]
        for Q in qs[1:]:
            perm=best_perm(Q, ref)
            aligned.append(Q[:,perm])
        meanQ=np.mean(np.stack(aligned,2), axis=2)
        meanQ = meanQ / meanQ.sum(axis=1, keepdims=True)
        # representative run = closest to mean
        d=[np.linalg.norm(A-meanQ) for A in aligned]
        rep_idx=int(np.argmin(d))
        # sort individuals by max membership
        meanQ_sorted, order = sort_inds(meanQ)
        rep_sorted = aligned[rep_idx][order,:]
        labels=[str(i+1) for i in range(meanQ.shape[0])]
        labels_sorted=[labels[i] for i in order]
        # save TSVs
        np.savetxt(PLOTS_DIR/f"Raccoon_STRUCTURE_meanQ_K{K}.tsv", meanQ, fmt="%.6f", delimiter="\t")
        np.savetxt(PLOTS_DIR/f"Raccoon_STRUCTURE_repQ_K{K}.tsv", aligned[rep_idx], fmt="%.6f", delimiter="\t")
        # plot
        barplot(meanQ_sorted, f"STRUCTURE — mean Q (K={K}, n={len(aligned)} runs)",
                labels_sorted, PLOTS_DIR/f"Raccoon_STRUCTURE_meanQ_K{K}.png", PLOTS_DIR/f"Raccoon_STRUCTURE_meanQ_K{K}.pdf")
        barplot(rep_sorted, f"STRUCTURE — representative run (K={K})\n{os.path.basename(keep[rep_idx])}",
                labels_sorted, PLOTS_DIR/f"Raccoon_STRUCTURE_repQ_K{K}.png", PLOTS_DIR/f"Raccoon_STRUCTURE_repQ_K{K}.pdf")
        print(f"[ok] K={K}: plots and tables written to {PLOTS_DIR}")

if __name__=="__main__":
    main()
PY
chmod 755 "$PLOTTER"

### -----------------------------
### RUN THE PLOTTER
### -----------------------------
echo "[info] Plotting STRUCTURE results..."
export ROOT_DIR  # so the Python script can pick it up
python3 "$PLOTTER"

echo
echo "[done] Outputs:"
echo "  - Harvester:  $HARVEST_DIR"
echo "  - Plots & TSV: $PLOTS_DIR"
