#!/usr/bin/env python3
"""
build_apkr_site.py — turn apkr_reactions.db into the single-file website
"The Asymmetric Pauson-Khand Database".

    python3 build_apkr_site.py apkr_reactions.db \
        --template apkr_site_template.html \
        --out apkr_database.html

What it does
------------
1. Reads every row of `all_data` (the merged table) plus the per-paper sheets,
   so each row can be traced to the spreadsheet it came from.
2. Pulls temperature, CO pressure, atmosphere and additives out of the free-text
   Notes column, where the extraction skill parks variables that have no column.
3. Featurizes each unique substrate and runs PCA (+ k-means) over the standardized
   descriptor matrix. Those coordinates are the chemical-space map on the front page.
4. Writes the data as JSON and, if a template is given, inlines it so the result
   is one HTML file that needs no web server.

The chemical-space map
---------------------
The map reproduces chemical_space_investigation.ipynb:

    Morgan fingerprint (radius 2, 1024 bits)
      -> PCA(n_components=80, random_state=42)
      -> TSNE(perplexity=20, max_iter=3000, learning_rate=150,
              early_exaggeration=15, init='pca', random_state=42)
      -> coordinates centred on the mean
      -> HDBSCAN(min_cluster_size=5, min_samples=5) on those 2D coordinates

PCA, t-SNE and HDBSCAN run here through scikit-learn with the notebook's own
hyperparameters. The fingerprint needs RDKit: if RDKit is importable the map is
bit-for-bit the notebook's, and if it is not, a surrogate circular fingerprint
from the bundled smiles_lite is used instead and every page says so.

Three ways to get the map, in order of fidelity:
  1. --coords cluster_MM.DD.YYYY.csv  reads the notebook's own output
     (columns: smiles, cluster, tsne_1, tsne_2) and uses those coordinates and
     cluster labels verbatim. No recomputation, nothing to reconcile.
  2. run this script where RDKit is installed - same pipeline, same seeds.
  3. no RDKit and no CSV - surrogate fingerprint, honest labelling.
"""

import argparse
import json
import math
import os
import re
import sqlite3
import sys
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Descriptors shown on substrate pages. They are not what the map is built
# from - the map uses the notebook's fingerprint pipeline - but they make a
# point on the map readable ("this is the big aromatic one").
FEATURES = [
    'heavy_atoms', 'mw', 'rings', 'aromatic_rings', 'aliphatic_rings',
    'rotatable_bonds', 'frac_csp3', 'frac_aromatic', 'heteroatoms',
    'n_N', 'n_O', 'n_Si', 'n_halogen', 'n_triple_bonds',
    'n_cc_double_bonds', 'n_carbonyl', 'defined_stereocentres',
]

# notebook hyperparameters, kept in one place
PCA_COMPONENTS = 80
PCA_SEED = 42
TSNE_KW = dict(n_components=2, perplexity=20, max_iter=3000, learning_rate=150,
               early_exaggeration=15, init='pca', random_state=42)
HDBSCAN_KW = dict(min_cluster_size=5, min_samples=5)
MORGAN_KW = dict(radius=2, fpSize=1024)


# ---------------------------------------------------------------- descriptors


def rdkit_features(smiles_list):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDLogger.DisableLog('rdApp.*')
    out = {}
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            out[smi] = None
            continue
        ri = m.GetRingInfo()
        arom_rings = sum(
            1 for r in ri.AtomRings()
            if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r))
        halo = sum(1 for a in m.GetAtoms() if a.GetSymbol() in ('F', 'Cl', 'Br', 'I'))
        triple = sum(1 for b in m.GetBonds()
                     if b.GetBondType() == Chem.BondType.TRIPLE)
        cc_double = sum(1 for b in m.GetBonds()
                        if b.GetBondType() == Chem.BondType.DOUBLE
                        and b.GetBeginAtom().GetSymbol() == 'C'
                        and b.GetEndAtom().GetSymbol() == 'C'
                        and not b.IsInRing())
        carbonyl = len(m.GetSubstructMatches(Chem.MolFromSmarts('[CX3]=[OX1]')))
        arom_atoms = sum(1 for a in m.GetAtoms() if a.GetIsAromatic())
        heavy = m.GetNumHeavyAtoms()
        out[smi] = {
            'heavy_atoms': heavy,
            'mw': round(Descriptors.MolWt(m), 3),
            'formula': rdMolDescriptors.CalcMolFormula(m),
            'rings': ri.NumRings(),
            'aromatic_rings': arom_rings,
            'aliphatic_rings': ri.NumRings() - arom_rings,
            'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(m),
            'frac_csp3': round(rdMolDescriptors.CalcFractionCSP3(m), 4),
            'frac_aromatic': round(arom_atoms / heavy, 4) if heavy else 0.0,
            'heteroatoms': rdMolDescriptors.CalcNumHeteroatoms(m),
            'n_N': sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'N'),
            'n_O': sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'O'),
            'n_S': sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'S'),
            'n_Si': sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'Si'),
            'n_halogen': halo,
            'n_triple_bonds': triple,
            'n_cc_double_bonds': cc_double,
            'n_carbonyl': carbonyl,
            'defined_stereocentres': len(Chem.FindMolChiralCenters(
                m, includeUnassigned=False, useLegacyImplementation=False)),
            'tpsa': round(Descriptors.TPSA(m), 2),
            'clogp': round(Descriptors.MolLogP(m), 3),
        }
    return out


def lite_features(smiles_list):
    import smiles_lite
    out = {}
    for smi in smiles_list:
        try:
            out[smi] = smiles_lite.descriptors(smi)
        except Exception as exc:                      # noqa: BLE001
            print(f'  ! could not parse {smi}: {exc}', file=sys.stderr)
            out[smi] = None
    return out


def have_rdkit():
    try:
        import rdkit  # noqa: F401
        return True
    except ImportError:
        return False


def featurize(smiles_list):
    if have_rdkit():
        return rdkit_features(smiles_list), 'rdkit'
    return lite_features(smiles_list), 'smiles_lite'


# ------------------------------------------------- the notebook's chemical map


def morgan_matrix(smiles_list):
    """Morgan fingerprints as an (n, 1024) float matrix. Returns (X, kept, how)."""
    if have_rdkit():
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
        RDLogger.DisableLog('rdApp.*')
        pairs = [(s, Chem.MolFromSmiles(s)) for s in smiles_list]
        pairs = [(s, m) for s, m in pairs if m is not None]
        gen = rdFingerprintGenerator.GetMorganGenerator(**MORGAN_KW)
        X = [gen.GetFingerprintAsNumPy(m).astype(float).tolist() for _, m in pairs]
        return X, [s for s, _ in pairs], 'rdkit Morgan r=2, 1024 bits'
    import smiles_lite
    X, kept = [], []
    for s in smiles_list:
        try:
            X.append([float(b) for b in smiles_lite.morgan_like_bits(
                s, radius=MORGAN_KW['radius'], nbits=MORGAN_KW['fpSize'])])
            kept.append(s)
        except Exception as exc:                      # noqa: BLE001
            print(f'  ! fingerprint failed for {s}: {exc}', file=sys.stderr)
    return X, kept, 'surrogate circular fingerprint r=2, 1024 bits (no RDKit)'


def chemical_space(smiles_list):
    """Morgan -> PCA(80) -> t-SNE -> HDBSCAN, exactly as in the notebook."""
    import numpy as np
    from sklearn.cluster import HDBSCAN
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE

    X, kept, how = morgan_matrix(smiles_list)
    X = np.asarray(X, dtype=np.float64)
    n = len(kept)
    ncomp = min(PCA_COMPONENTS, n, X.shape[1])
    pca = PCA(n_components=ncomp, random_state=PCA_SEED)
    X_pca = pca.fit_transform(X)
    kw = dict(TSNE_KW)
    kw['perplexity'] = min(kw['perplexity'], max(5, (n - 1) // 3))
    X_tsne = TSNE(**kw).fit_transform(X_pca)
    X_tsne = X_tsne - X_tsne.mean(axis=0)
    try:
        labels = HDBSCAN(copy=True, **HDBSCAN_KW).fit_predict(X_tsne)
    except TypeError:                      # older scikit-learn has no copy kwarg
        labels = HDBSCAN(**HDBSCAN_KW).fit_predict(X_tsne)
    print(f'  {how}: {n} substrates, PCA({ncomp}) keeps '
          f'{pca.explained_variance_ratio_.sum():.0%}, '
          f'{len(set(labels) - {-1})} HDBSCAN clusters, '
          f'{int((labels == -1).sum())} noise points')
    return {
        'smiles': kept,
        'tsne': [[round(float(a), 4), round(float(b), 4)] for a, b in X_tsne],
        'pc': [[round(float(r[0]), 4), round(float(r[1]), 4)] for r in X_pca],
        'labels': [int(v) for v in labels],
        'how': how,
        'pca_components': int(ncomp),
        'pca_explained': round(float(pca.explained_variance_ratio_.sum()), 4),
        'pc12_explained': [round(float(v), 4)
                           for v in pca.explained_variance_ratio_[:2]],
        'source': 'recomputed',
    }


def chemical_space_from_csv(path):
    """Read the notebook's own cluster_MM.DD.YYYY.csv instead of recomputing."""
    import csv
    rows = list(csv.DictReader(open(path, newline='')))
    if not rows:
        sys.exit(f'{path} is empty')
    keys = {k.lower().strip(): k for k in rows[0]}

    def pick(*names):
        for nm in names:
            if nm in keys:
                return keys[nm]
        return None
    c_smi = pick('smiles', 'smile', 'canonical_smiles')
    c_x = pick('tsne_1', 'tsne1', 'tnse_1', 'x')
    c_y = pick('tsne_2', 'tnse_2', 'tsne2', 'y')
    c_lab = pick('cluster', 'label', 'labels')
    if not (c_smi and c_x and c_y):
        sys.exit(f'{path} needs smiles, tsne_1 and tsne_2 columns; found {list(rows[0])}')
    out = {'smiles': [], 'tsne': [], 'pc': [], 'labels': [],
           'how': 'coordinates read from ' + os.path.basename(path),
           'pca_components': PCA_COMPONENTS, 'pca_explained': None,
           'pc12_explained': None, 'source': 'csv'}
    for r in rows:
        out['smiles'].append(r[c_smi].strip())
        out['tsne'].append([round(float(r[c_x]), 4), round(float(r[c_y]), 4)])
        out['labels'].append(int(float(r[c_lab])) if c_lab and r[c_lab] != '' else -1)
    xs = [p[0] for p in out['tsne']]
    ys = [p[1] for p in out['tsne']]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    out['tsne'] = [[round(p[0] - mx, 4), round(p[1] - my, 4)] for p in out['tsne']]
    out['pc'] = out['tsne']
    print(f"  {out['how']}: {len(out['smiles'])} substrates, "
          f"{len(set(out['labels']) - {-1})} clusters, "
          f"{out['labels'].count(-1)} noise points")
    return out


# ------------------------------------------------------------------ text mining

TEMP_PATTERNS = [
    # "100 +/- 3 C oil bath" -> 100, not 3
    r'(-?\d+(?:\.\d+)?)\s*(?:\+/-|\+-|±)\s*\d+(?:\.\d+)?\s*(?:°|º)?\s*C\b',
    r'[Tt]emp(?:erature)?[^;.\n]*?(-?\d+(?:\.\d+)?)\s*(?:°|º)?\s*C\b',
    r'(-?\d+(?:\.\d+)?)\s*(?:°|º)\s*C\b',
    r'\bat\s+(-?\d+(?:\.\d+)?)\s*C\b',
    r'(-?\d+(?:\.\d+)?)\s*C\b',
]


def parse_temp(notes):
    if not notes:
        return None
    for pat in TEMP_PATTERNS:
        m = re.search(pat, notes)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
    if re.search(r'\b(rt|r\.t\.|room temp)', notes, re.I):
        return 25.0
    return None


def parse_field(notes, label):
    if not notes:
        return None
    m = re.search(label + r'\s*([^;.\n]+)', notes, re.I)
    return m.group(1).strip() if m else None


def first_number(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    m = re.search(r'-?\d+(?:\.\d+)?', str(value))
    return float(m.group()) if m else None


def norm(v):
    """Stringify a cell so that 57, 57.0 and '57.0' all compare equal."""
    if v is None:
        return ''
    if isinstance(v, (int, float)):
        return str(int(v)) if float(v) == int(v) else str(float(v))
    s = str(v).strip()
    try:
        f = float(s)
    except ValueError:
        return s
    return str(int(f)) if f == int(f) else str(f)


JOURNALS = {
    'JACS': 'J. Am. Chem. Soc.', 'JOC': 'J. Org. Chem.',
    'Organometallics': 'Organometallics', 'Tetrahedron': 'Tetrahedron',
    'Adv-Synth-Catal': 'Adv. Synth. Catal.', 'AdvSynthCatal': 'Adv. Synth. Catal.',
    'Adv_Synth_Catal': 'Adv. Synth. Catal.', 'ChemCommun': 'Chem. Commun.',
    'ChemAsianJ': 'Chem. Asian J.', 'ACS-Catal': 'ACS Catal.',
    'Tetrahedron_Lett': 'Tetrahedron Lett.', 'OCF': 'Org. Chem. Front.',
    'AngewChemIntEd': 'Angew. Chem. Int. Ed.',
}


def sheet_label(sheet):
    """'Kim_JOC_2008_08.07.2026.xlsx' -> ('Kim', 'J. Org. Chem.', '2008')"""
    stem = re.sub(r'\.xlsx$', '', sheet)
    parts = stem.split('_')
    author = parts[0]
    year = None
    journal = None
    for p in parts[1:]:
        if re.fullmatch(r'(19|20)\d{2}', p):
            year = p
            break
        journal = p if journal is None else journal + '_' + p
    if year is None:
        m = re.search(r'(19|20)\d{2}', stem)
        year = m.group() if m else ''
    jr = JOURNALS.get(journal or '', (journal or '').replace('-', ' ').replace('_', ' '))
    return author.replace('-', ' '), jr, year


# ----------------------------------------------------------------------- build


def build(db_path, coords_csv=None):
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    cur = con.cursor()

    db_cols = [r[1] for r in cur.execute('PRAGMA table_info(all_data)')]
    base_cols = db_cols[:27]

    sheets = [r[0] for r in cur.execute('SELECT sheet_name FROM added_sheets')]

    # --- row -> source sheet attribution, by matching the 27 shared columns
    sel = ','.join(f'"{c}"' for c in base_cols)
    sheet_index = {}
    for sh in sheets:
        have = {r[1] for r in cur.execute(f'PRAGMA table_info("{sh}")')}
        if not set(base_cols) <= have:
            continue
        for row in cur.execute(f'SELECT {sel} FROM "{sh}"'):
            sheet_index.setdefault(tuple(norm(v) for v in row), []).append(sh)

    rows = [dict(r) for r in cur.execute('SELECT * FROM all_data')]

    # --- papers
    doi_sheets = {}
    for sh in sheets:
        try:
            for (doi,) in cur.execute(f'SELECT DISTINCT DOI FROM "{sh}"'):
                doi_sheets.setdefault(norm(doi), set()).add(sh)
        except sqlite3.Error:
            pass

    papers = []
    paper_idx = {}
    for doi in sorted({norm(r['DOI']) for r in rows}):
        shs = sorted(doi_sheets.get(doi, []))
        author = journal = year = ''
        for sh in shs:
            a, j, y = sheet_label(sh)
            if a and not re.fullmatch(r'(19|20)\d{2}', a) and (not author or len(a) < len(author)):
                if j:
                    author, journal, year = a, j, y
        if not author and shs:
            author, journal, year = sheet_label(shs[0])
        paper_idx[doi] = len(papers)
        papers.append({
            'doi': doi,
            'author': author,
            'journal': journal,
            'year': year,
            'label': f'{author} {year}'.strip() if author else doi,
            'citation': ', '.join(x for x in (author, journal, year) if x) or doi,
            'sheets': shs,
            'n_rows': 0,
        })

    # --- substrates and ligands
    sub_smiles, lig_smiles = [], []
    for r in rows:
        s = norm(r['canonical_sub_smiles']) or norm(r['Substrate SMILES'])
        if s and s not in sub_smiles:
            sub_smiles.append(s)
        l = norm(r['ligand_smiles']) or norm(r['Ligand SMILES'])
        if l and l not in lig_smiles:
            lig_smiles.append(l)

    feats, engine = featurize(sub_smiles)
    lig_feats, _ = featurize(lig_smiles)

    space = (chemical_space_from_csv(coords_csv) if coords_csv
             else chemical_space(sub_smiles))
    space_index = {s: j for j, s in enumerate(space['smiles'])}

    sub_index = {s: i for i, s in enumerate(sub_smiles)}
    substrates = []
    unplaced = []
    for i, s in enumerate(sub_smiles):
        d = feats.get(s) or {}
        j = space_index.get(s)
        substrates.append({
            'id': i,
            'smiles': s,
            'formula': d.get('formula', ''),
            'desc': {f: d.get(f) for f in FEATURES if f in d},
            'tsne': space['tsne'][j] if j is not None else None,
            'pc': space['pc'][j] if j is not None else None,
            'cluster': space['labels'][j] if j is not None else None,
            'rows': [],
            'cpd': [],
            'papers': [],
        })
        if j is None:
            unplaced.append(s)
    if unplaced:
        print(f'  ! {len(unplaced)} substrate(s) have no coordinates and will not '
              f'appear on the map:', file=sys.stderr)
        for s in unplaced:
            print('    ', s, file=sys.stderr)

    lig_index = {s: i for i, s in enumerate(lig_smiles)}
    ligands = []
    for i, s in enumerate(lig_smiles):
        d = lig_feats.get(s) or {}
        ligands.append({
            'id': i, 'smiles': s, 'formula': d.get('formula', ''),
            'names': [], 'rows': [], 'papers': [], 'cpd': [],
            'mw': d.get('mw'), 'heavy_atoms': d.get('heavy_atoms'),
        })

    # --- rows
    out_rows, derived = [], []
    for rid, r in enumerate(rows):
        doi = norm(r['DOI'])
        s = norm(r['canonical_sub_smiles']) or norm(r['Substrate SMILES'])
        l = norm(r['ligand_smiles']) or norm(r['Ligand SMILES'])
        si = sub_index.get(s, -1)
        li = lig_index.get(l, -1)
        pi = paper_idx.get(doi, -1)
        notes = norm(r['Notes'])
        src = sheet_index.get(tuple(norm(r[c]) for c in base_cols), [])

        out_rows.append([norm(r[c]) for c in db_cols])
        derived.append({
            'sub': si, 'lig': li, 'paper': pi,
            'temp': parse_temp(notes),
            'co': parse_field(notes, r'CO pressure'),
            'atm': parse_field(notes, r'atmosphere'),
            'add': parse_field(notes, r'additives?'),
            'yield': first_number(r['% yield (product)']),
            'ee': first_number(r['% ee (product)']),
            'src': src,
        })
        if pi >= 0:
            papers[pi]['n_rows'] += 1
        if si >= 0:
            substrates[si]['rows'].append(rid)
            cpd = norm(r['Substrate cpd #'])
            if cpd and cpd not in substrates[si]['cpd']:
                substrates[si]['cpd'].append(cpd)
            if pi >= 0 and pi not in substrates[si]['papers']:
                substrates[si]['papers'].append(pi)
        if li >= 0:
            ligands[li]['rows'].append(rid)
            name = norm(r['ligand_preferred_name']) or norm(r['Ligand name'])
            if name and name not in ligands[li]['names']:
                ligands[li]['names'].append(name)
            cpd = norm(r['Ligand cpd #'])
            if cpd and cpd not in ligands[li]['cpd']:
                ligands[li]['cpd'].append(cpd)
            if pi >= 0 and pi not in ligands[li]['papers']:
                ligands[li]['papers'].append(pi)

    # --- repeat-extraction detection: same paper + same reported numbers, but
    #     the rows came from two different spreadsheets (a re-parse of one paper)
    groups = {}
    for rid, d in enumerate(derived):
        r = rows[rid]
        key = (d['paper'], norm(r['Substrate cpd #']), norm(r['Product cpd #']),
               norm(r['% yield (product)']), norm(r['% ee (product)']),
               norm(r['ligand_preferred_name']) or norm(r['Ligand name']))
        groups.setdefault(key, []).append(rid)
    dup_group = [None] * len(derived)
    gid = 0
    for key, ids in groups.items():
        if len(ids) < 2:
            continue
        srcs = {tuple(derived[i]['src']) for i in ids}
        if len(srcs) < 2:
            continue                      # genuinely repeated entry inside one sheet
        for i in ids:
            dup_group[i] = gid
        gid += 1
    for rid, d in enumerate(derived):
        d['dup'] = dup_group[rid]

    for s in substrates:
        ees = [derived[i]['ee'] for i in s['rows'] if derived[i]['ee'] is not None]
        ys = [derived[i]['yield'] for i in s['rows'] if derived[i]['yield'] is not None]
        s['best_ee'] = max(ees) if ees else None
        s['best_yield'] = max(ys) if ys else None
        s['n_rows'] = len(s['rows'])
    for l in ligands:
        l['n_rows'] = len(l['rows'])

    data = {
        'meta': {
            'title': 'The Asymmetric Pauson-Khand Database',
            'generated': date.today().strftime('%d %B %Y'),
            'source_db': os.path.basename(db_path),
            'n_rows': len(out_rows),
            'n_papers': len(papers),
            'n_substrates': len(substrates),
            'n_ligands': len(ligands),
            'n_repeat_rows': sum(1 for d in derived if d['dup'] is not None),
            'featurizer': engine,
            'features': FEATURES,
            'space': {
                'how': space['how'],
                'source': space['source'],
                'exact': space['source'] == 'csv' or engine == 'rdkit',
                'n_placed': len(space['smiles']),
                'n_unplaced': len(unplaced),
                'pca_components': space['pca_components'],
                'pca_explained': space['pca_explained'],
                'pc12_explained': space['pc12_explained'],
                'clusters': sorted(set(space['labels']) - {-1}),
                'n_noise': space['labels'].count(-1),
                'tsne_params': {k: v for k, v in TSNE_KW.items()},
                'hdbscan_params': HDBSCAN_KW,
                'morgan_params': MORGAN_KW,
            },
        },
        'columns': db_cols,
        'rows': out_rows,
        'derived': derived,
        'papers': papers,
        'substrates': substrates,
        'ligands': ligands,
    }
    con.close()
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('db')
    ap.add_argument('--template', default=None)
    ap.add_argument('--out', default='apkr_database.html')
    ap.add_argument('--json', default=None)
    ap.add_argument('--coords', default=None,
                    help="the notebook's cluster_MM.DD.YYYY.csv; its t-SNE "
                         'coordinates and HDBSCAN labels are used verbatim')
    args = ap.parse_args()

    data = build(args.db, args.coords)
    blob = json.dumps(data, separators=(',', ':'), ensure_ascii=False)
    # '<' only ever appears inside JSON strings here, so escaping it wholesale is
    # safe and stops a stray '</script' in a Notes cell from ending the data block
    blob = blob.replace('<', '\\u003c').replace('\u2028', '\\u2028').replace('\u2029', '\\u2029')
    m = data['meta']
    print(f"{m['n_rows']} reactions | {m['n_papers']} papers | "
          f"{m['n_substrates']} substrates | {m['n_ligands']} ligands")
    if not m['space']['exact']:
        print('NOTE: built without RDKit, so map coordinates are a surrogate. '
              'Rerun where RDKit lives, or pass --coords cluster_MM.DD.YYYY.csv.')

    if args.json:
        with open(args.json, 'w') as fh:
            fh.write(blob)
        print('wrote', args.json)

    if args.template:
        with open(args.template) as fh:
            html = fh.read()
        token = '"__APKR_DATA__"'
        if token not in html:
            sys.exit(f'template {args.template} has no {token} placeholder')
        html = html.replace(token, blob)
        with open(args.out, 'w') as fh:
            fh.write(html)
        print('wrote', args.out, f'({len(html) / 1e6:.2f} MB)')


if __name__ == '__main__':
    main()
