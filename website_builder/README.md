# The Asymmetric Pauson–Khand Database — website

A single HTML file that reads nothing and needs no server: open
`apkr_database.html` in a browser and the whole database is there. Built
22 September 2026 from `apkr_reactions.db` — 716 reactions, 136 unique
substrates, 39 unique ligands, 16 papers.

## Files

| file | what it is |
|---|---|
| `apkr_database.html` | the site, data included (1.1 MB). Double-click it. |
| `build_apkr_site.py` | rebuilds the site from the database |
| `apkr_site_template.html` | the app itself; the builder injects data into it |
| `smiles_lite.py` | dependency-free SMILES reader, used only when RDKit is missing |

## Rebuilding after you add a paper

```bash
python3 build_apkr_site.py apkr_reactions.db \
    --template apkr_site_template.html \
    --out apkr_database.html
```

Run it from a place where RDKit and scikit-learn are importable — the same
environment as `chemical_space_investigation.ipynb`. It reports what it did:

```
rdkit Morgan r=2, 1024 bits: 136 substrates, PCA(80) keeps 99%,
  8 HDBSCAN clusters, 2 noise points
716 reactions | 16 papers | 136 substrates | 39 ligands
```

If it prints a note about a surrogate fingerprint, RDKit was not importable and
the map is provisional — the site says so on the map page too.

### Using the notebook's own coordinates instead

```bash
python3 build_apkr_site.py apkr_reactions.db --coords cluster_08.13.2026.csv \
    --template apkr_site_template.html --out apkr_database.html
```

The t-SNE coordinates and HDBSCAN labels in that CSV are then used verbatim, so
the map matches your notebook figure and the per-cluster grid images exactly,
with nothing recomputed. It reads either `tsne_2` or the notebook's `tnse_2`
spelling.

## What the map is

The pipeline from the notebook, unchanged in hyperparameters and seeds:
Morgan fingerprint (radius 2, 1024 bits) → `PCA(n_components=80,
random_state=42)` → `TSNE(perplexity=20, max_iter=3000, learning_rate=150,
early_exaggeration=15, init='pca', random_state=42)` → centre on the mean →
`HDBSCAN(min_cluster_size=5, min_samples=5)`.

Cluster numbers are HDBSCAN's own labels, so they line up with your
`cluster_*.csv`. Points HDBSCAN called noise (−1) are drawn grey and labelled
unclustered rather than dropped.

One point per reaction by default, so a substrate used 39 times shows 39 points
spread into a small rosette around its coordinates — that is what makes hover
useful, since each point carries its own yield, ee and conditions. Switch to one
point per substrate in the sidebar when you want the map to read as chemical
space alone. The axis selector also plots PC1 vs PC2 of the same fingerprint
matrix, worth a look because distances between t-SNE clusters mean nothing while
PCA distances do.

## Structure search

Substrate similarity and the ligand check use a small SMILES reader built into
the page, not RDKit: your typed query and the stored entries are fingerprinted
by the same code, with no server round-trip. It hashes circular atom
environments to radius 2 (ECFP-like, not RDKit-identical), so treat it as "is
this skeleton already here", not as a number to quote. **It ignores
stereochemistry** — a hit tells you the skeleton is present, not that the
enantiomer or the P-stereogenic variant is. A BINAP SMILES written in a
completely different atom order still matches the stored BINAP at Tanimoto
1.000.

## Two things in the data worth knowing

**150 rows are the same published entry extracted twice.** The Jeong
*Chem. Commun.* 2004 sheets from 08.16 and 08.19 hold identical numbers with
differently written SMILES and notes, and `Qi_ACS-Catal_2024` overlaps heavily
with `Qi-2024-redo-rdkit`. Rows are flagged `repeat` wherever they appear, each
entry page links to its twin, and the map has a switch to hide them. Nothing was
deleted. If you'd rather resolve it properly, the two Jeong sheets are exact
duplicates and one can be dropped from `all_data`; the Qi pair needs a decision
about which parse is canonical.

**Temperature, CO pressure, atmosphere and additives are parsed from the Notes
text**, since the template has no columns for them. Temperature resolves for 669
of 716 rows. The parser prefers the real value in phrasings like
`100 +/- 3 C oil bath`, but it is still regex over prose, so the full note sits
on every entry page. If these become columns in the template, the parsing drops
out of `build_apkr_site.py` and the fields come straight from the database.

## Adding a field to the entry pages

Nothing to edit: entry pages iterate over whatever columns `all_data` has. Add a
column to the template and the spreadsheets, rebuild, and it appears. The four
parsed condition fields and the hover card are the only places that name
specific columns.

## Known limits

- Structure drawings come from the smiles-drawer library, loaded from a CDN on
  first open — the only thing here that touches the network. Offline you get
  formulas and SMILES text instead. To make it fully offline, drop
  `smiles-drawer.min.js` next to the HTML; the page looks for a local copy first.
- The site was verified headlessly (every route, all 716 entry pages, both
  search tools, the CSV export) but not in a real browser, so give the drawings
  and the hover card a look on your machine.
- 39 ligand SMILES cover 694 of 716 rows; the remaining 22 record a ligand name
  with no SMILES (ligand-free runs, or "not specified in this paper"). They are
  in the tables but not in the ligand check.
