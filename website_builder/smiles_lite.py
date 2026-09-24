"""
smiles_lite.py — a small, dependency-free SMILES reader used to featurize the
APKR substrate set inside a sandbox that has no RDKit installed.

It parses the subset of SMILES that appears in this database (organic subset,
bracket atoms, ring-closure digits and %nn, branches, bond symbols, aromatic
lowercase atoms, charges, isotopes, stereo markers) into an atom/bond graph and
derives interpretable descriptors from that graph.

It is deliberately NOT a chemistry toolkit: no aromaticity perception, no CIP,
no ring perception beyond the cycle rank / smallest-set heuristics documented
below. Every value it produces is a graph or composition count, which is what
the PCA in build_apkr_site.py consumes. If RDKit is importable, that script
prefers RDKit and this module is bypassed.
"""

import re

# average atomic masses (enough elements for this dataset)
MASS = {
    'H': 1.008, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
    'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Br': 79.904, 'I': 126.904,
    'Se': 78.971, 'Na': 22.990, 'K': 39.098, 'Li': 6.94, 'Zn': 65.38, 'Ag': 107.868,
    'Rh': 102.906, 'Ir': 192.217, 'Co': 58.933, 'Pd': 106.42, 'Ti': 47.867,
    'Mo': 95.95, 'W': 183.84, 'Fe': 55.845, 'Cu': 63.546, 'Al': 26.982,
}

# allowed valences, smallest first; implicit H fills to the first valence >= bond sum
VALENCE = {
    'B': [3], 'C': [4], 'N': [3, 5], 'O': [2], 'F': [1], 'Si': [4], 'P': [3, 5],
    'S': [2, 4, 6], 'Cl': [1], 'Br': [1], 'I': [1], 'Se': [2, 4, 6],
}

ORGANIC = ['Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I',
           'c', 'n', 'o', 'p', 's', 'b']

HALOGENS = {'F', 'Cl', 'Br', 'I'}

_BRACKET = re.compile(
    r'\[(?P<iso>\d+)?(?P<sym>[A-Za-z][a-z]?|\*)(?P<chiral>@{1,2}(?:TH[12]|AL[12]|SP[123]|TB\d{1,2}|OH\d{1,2})?)?'
    r'(?P<hyd>H\d?)?(?P<chg>(?:\+{1,3}|-{1,3}|\+\d|-\d))?(?::\d+)?\]'
)


class Atom:
    __slots__ = ('sym', 'aromatic', 'charge', 'explicit_h', 'chiral', 'nbrs',
                 'bond_sum', 'in_ring', 'idx')

    def __init__(self, sym, aromatic, idx):
        self.sym = sym
        self.aromatic = aromatic
        self.charge = 0
        self.explicit_h = None
        self.chiral = False
        self.nbrs = []          # (atom_idx, order)
        self.bond_sum = 0.0
        self.in_ring = False
        self.idx = idx


class ParseError(ValueError):
    pass


def parse(smiles):
    """Return (atoms, bonds) where bonds is a list of (i, j, order, in_ring)."""
    atoms = []
    bonds = []
    stack = []
    ring = {}                  # ring-closure label -> (atom idx, pending order)
    prev = None
    pending_order = None
    i = 0
    n = len(smiles)
    s = smiles

    def add_bond(a, b, order):
        bonds.append([a, b, order, False])
        atoms[a].nbrs.append((b, order))
        atoms[b].nbrs.append((a, order))
        atoms[a].bond_sum += order
        atoms[b].bond_sum += order

    while i < n:
        ch = s[i]
        if ch == '[':
            m = _BRACKET.match(s, i)
            if not m:
                raise ParseError(f'bad bracket atom at {i} in {s}')
            sym = m.group('sym')
            arom = sym[:1].islower() and sym not in ('Cl', 'Br')
            a = Atom(sym.capitalize() if arom else sym, arom, len(atoms))
            if m.group('chg'):
                c = m.group('chg')
                if c[0] == '+':
                    a.charge = int(c[1:]) if c[1:].isdigit() else len(c)
                else:
                    a.charge = -(int(c[1:]) if c[1:].isdigit() else len(c))
            h = m.group('hyd')
            a.explicit_h = 0 if h is None else (1 if h == 'H' else int(h[1:]))
            a.chiral = bool(m.group('chiral'))
            atoms.append(a)
            if prev is not None:
                add_bond(prev, a.idx, pending_order if pending_order else 1)
            prev = a.idx
            pending_order = None
            i = m.end()
            continue
        matched = None
        for sym in ORGANIC:
            if s.startswith(sym, i):
                matched = sym
                break
        if matched:
            arom = matched.islower()
            a = Atom(matched.upper() if arom else matched, arom, len(atoms))
            atoms.append(a)
            if prev is not None:
                add_bond(prev, a.idx, pending_order if pending_order else (1.5 if (arom and atoms[prev].aromatic) else 1))
            prev = a.idx
            pending_order = None
            i += len(matched)
            continue
        if ch in '-=#$:/\\':
            pending_order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5,
                             '/': 1, '\\': 1}[ch]
            i += 1
            continue
        if ch == '(':
            stack.append(prev)
            i += 1
            continue
        if ch == ')':
            prev = stack.pop()
            i += 1
            continue
        if ch == '%':
            label = s[i + 1:i + 3]
            i += 3
            _close_ring(atoms, ring, label, prev, pending_order, add_bond)
            pending_order = None
            continue
        if ch.isdigit():
            _close_ring(atoms, ring, ch, prev, pending_order, add_bond)
            pending_order = None
            i += 1
            continue
        if ch == '.':
            prev = None
            pending_order = None
            i += 1
            continue
        if ch in '*~':
            i += 1
            continue
        raise ParseError(f'unhandled char {ch!r} at {i} in {s}')

    if ring:
        raise ParseError(f'unclosed ring bond(s) {sorted(ring)} in {s}')
    _mark_ring_bonds(atoms, bonds)
    return atoms, bonds


def _close_ring(atoms, ring, label, prev, pending_order, add_bond):
    if prev is None:
        raise ParseError('ring closure with no preceding atom')
    if label in ring:
        other, order = ring.pop(label)
        o = pending_order or order
        if o is None:
            o = 1.5 if (atoms[prev].aromatic and atoms[other].aromatic) else 1
        add_bond(other, prev, o)
    else:
        ring[label] = (prev, pending_order)


def _mark_ring_bonds(atoms, bonds):
    """Flag every bond that lies on a cycle (bridge detection, iterative DFS)."""
    n = len(atoms)
    adj = [[] for _ in range(n)]
    for bi, (a, b, o, _) in enumerate(bonds):
        adj[a].append((b, bi))
        adj[b].append((a, bi))
    disc = [0] * n
    low = [0] * n
    seen = [False] * n
    timer = [1]
    bridges = set()
    for root in range(n):
        if seen[root]:
            continue
        stack = [(root, -1, iter(adj[root]))]
        seen[root] = True
        disc[root] = low[root] = timer[0]
        timer[0] += 1
        while stack:
            v, pbond, it = stack[-1]
            advanced = False
            for w, bi in it:
                if bi == pbond:
                    continue
                if not seen[w]:
                    seen[w] = True
                    disc[w] = low[w] = timer[0]
                    timer[0] += 1
                    stack.append((w, bi, iter(adj[w])))
                    advanced = True
                    break
                low[v] = min(low[v], disc[w])
            if not advanced:
                stack.pop()
                if stack:
                    u, ub, _ = stack[-1]
                    low[u] = min(low[u], low[v])
                    if low[v] > disc[u]:
                        bridges.add(pbond)
    for bi, bond in enumerate(bonds):
        if bi not in bridges:
            bond[3] = True
            atoms[bond[0]].in_ring = True
            atoms[bond[1]].in_ring = True


def implicit_h(a):
    if a.explicit_h is not None:
        return a.explicit_h
    vs = VALENCE.get(a.sym)
    if not vs:
        return 0
    need = a.bond_sum
    if a.aromatic:
        # an aromatic ring atom written lowercase: treat the ring bonds as 1.5
        need = max(need, 3.0 if a.sym in ('C', 'N') else 2.0)
    target = None
    for v in vs:
        if v + a.charge >= need - 1e-6:
            target = v + a.charge
            break
    if target is None:
        target = vs[-1] + a.charge
    return max(0, int(round(target - need)))


def _rotatable(atoms, bonds):
    """Acyclic single bonds between two non-terminal heavy atoms, minus amide C-N."""
    rot = 0
    for a, b, o, in_ring in bonds:
        if in_ring or o != 1:
            continue
        A, B = atoms[a], atoms[b]
        if len(A.nbrs) < 2 or len(B.nbrs) < 2:
            continue
        # RDKit's strict definition: a single bond next to a triple bond is
        # part of a linear fragment, so it does not count as rotatable
        if any(order == 3 for _, order in A.nbrs) or any(order == 3 for _, order in B.nbrs):
            continue
        # RDKit's strict definition also drops the C(=O)-N of an amide and the
        # C(=O)-O of an ester (delocalised, not freely rotating). The S-N of a
        # sulfonamide *is* counted, matching RDKit.
        pair = {A.sym, B.sym}
        if pair in ({'C', 'N'}, {'C', 'O'}):
            c = A if A.sym == 'C' else B
            if any(atoms[j].sym == 'O' and order == 2 for j, order in c.nbrs):
                continue
        rot += 1
    return rot


def _ring_count(atoms, bonds):
    """Cycle rank = bonds - atoms + connected components (SSSR size)."""
    n = len(atoms)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b, o, r in bonds:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = len({find(i) for i in range(n)}) if n else 0
    return len(bonds) - n + comp


def _aromatic_rings(atoms, bonds):
    """Count of 5/6-membered all-aromatic rings, found by walking ring bonds."""
    arom_bonds = [(a, b) for a, b, o, r in bonds
                  if r and atoms[a].aromatic and atoms[b].aromatic]
    if not arom_bonds:
        return 0
    adj = {}
    for a, b in arom_bonds:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    # cycle rank of the aromatic subgraph = number of aromatic rings
    nodes = set(adj)
    parent = {x: x for x in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in arom_bonds:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = len({find(x) for x in nodes})
    return len(arom_bonds) - len(nodes) + comp


def descriptors(smiles):
    atoms, bonds = parse(smiles)
    heavy = len(atoms)
    counts = {}
    for a in atoms:
        counts[a.sym] = counts.get(a.sym, 0) + 1
    nh = sum(implicit_h(a) for a in atoms)
    mw = sum(MASS.get(a.sym, 0.0) for a in atoms) + nh * MASS['H']

    carbons = [a for a in atoms if a.sym == 'C']
    sp3 = 0
    for a in carbons:
        if a.aromatic:
            continue
        if any(order > 1.4 for _, order in a.nbrs):
            continue
        sp3 += 1
    rings = _ring_count(atoms, bonds)
    arom_rings = _aromatic_rings(atoms, bonds)
    triple = sum(1 for a, b, o, r in bonds if o == 3)
    dbl_cc = sum(1 for a, b, o, r in bonds
                 if o == 2 and atoms[a].sym == 'C' and atoms[b].sym == 'C' and not r)
    carbonyl = sum(1 for a, b, o, r in bonds
                   if o == 2 and {atoms[a].sym, atoms[b].sym} == {'C', 'O'})
    halogen = sum(v for k, v in counts.items() if k in HALOGENS)
    hetero = sum(v for k, v in counts.items() if k not in ('C',))
    stereo = sum(1 for a in atoms if a.chiral)
    arom_atoms = sum(1 for a in atoms if a.aromatic)

    return {
        'heavy_atoms': heavy,
        'mw': round(mw, 3),
        'n_H': nh,
        'formula': _formula(counts, nh),
        'rings': rings,
        'aromatic_rings': arom_rings,
        'aliphatic_rings': max(0, rings - arom_rings),
        'rotatable_bonds': _rotatable(atoms, bonds),
        'frac_csp3': round(sp3 / len(carbons), 4) if carbons else 0.0,
        'frac_aromatic': round(arom_atoms / heavy, 4) if heavy else 0.0,
        'heteroatoms': hetero,
        'n_N': counts.get('N', 0),
        'n_O': counts.get('O', 0),
        'n_S': counts.get('S', 0),
        'n_Si': counts.get('Si', 0),
        'n_halogen': halogen,
        'n_triple_bonds': triple,
        'n_cc_double_bonds': dbl_cc,
        'n_carbonyl': carbonyl,
        'defined_stereocentres': stereo,
    }


def _fnv(s):
    h = 0x811c9dc5
    for ch in s:
        h ^= ord(ch)
        h = (h * 0x01000193) & 0xFFFFFFFF
    return h


def morgan_like_bits(smiles, radius=2, nbits=1024):
    """A circular (ECFP-style) bit fingerprint, folded to `nbits`.

    This mirrors CHEM.fingerprint() in the website exactly, so similarity
    computed in the browser agrees with anything computed here. It is NOT
    RDKit's Morgan fingerprint: the atom invariants and the hash differ, so
    bit indices are not comparable with RDKit's. Stereochemistry is ignored.
    """
    atoms, bonds = parse(smiles)
    n = len(atoms)
    bits = [0] * nbits
    cur = []
    for a in atoms:
        cur.append(_fnv('|'.join(str(x) for x in (
            a.sym, len(a.nbrs), implicit_h(a), a.charge,
            1 if a.aromatic else 0, 1 if a.in_ring else 0))))
    for h in cur:
        bits[h % nbits] = 1
    for r in range(1, radius + 1):
        nxt = []
        for i in range(n):
            env = sorted('%s:%s' % (o, cur[j]) for j, o in atoms[i].nbrs)
            nxt.append(_fnv('%d|%d|%s' % (r, cur[i], ','.join(env))))
        cur = nxt
        for h in cur:
            bits[h % nbits] = 1
    return bits


def _formula(counts, nh):
    order = ['C', 'H', 'N', 'O', 'S', 'Si', 'P', 'B', 'F', 'Cl', 'Br', 'I', 'Se']
    parts = []
    c = dict(counts)
    if nh:
        c['H'] = c.get('H', 0) + nh
    for sym in order:
        if c.get(sym):
            parts.append(sym + (str(c[sym]) if c[sym] > 1 else ''))
            c.pop(sym)
    for sym in sorted(c):
        if c[sym]:
            parts.append(sym + (str(c[sym]) if c[sym] > 1 else ''))
    return ''.join(parts)


if __name__ == '__main__':
    tests = {
        'C#CC(OCC=C)OCC=C': dict(heavy_atoms=11, rings=0, rotatable_bonds=6,
                                 frac_csp3=1 / 3, formula='C9H12O2'),
        'c1ccccc1': dict(heavy_atoms=6, rings=1, aromatic_rings=1, formula='C6H6'),
        'CC(=O)OC=C=C(C)CCCC#C[Si](C)(C)C': dict(heavy_atoms=17, rings=0,
                                                 formula='C14H22O2Si'),
        'C=C(C)CN(CC#CCl)S(=O)(=O)c1ccc(C)cc1': dict(heavy_atoms=19, rings=1,
                                                     aromatic_rings=1),
        'O=C1C(C)=C2CN(S(=O)(=O)c3ccc(C)cc3)[C@H](C=C)[C@@H]2C1': dict(rings=3,
                                                                       aromatic_rings=1,
                                                                       defined_stereocentres=2),
    }
    for smi, exp in tests.items():
        d = descriptors(smi)
        bad = {}
        for k, v in exp.items():
            got = d[k]
            if isinstance(v, str):
                if got != v:
                    bad[k] = (v, got)
            elif abs(got - v) > 1e-3:
                bad[k] = (v, got)
        print(('OK  ' if not bad else 'FAIL'), smi)
        if bad:
            print('     expected/got:', bad)
        print('     ', {k: d[k] for k in ('formula', 'mw', 'heavy_atoms', 'rings',
                                          'aromatic_rings', 'rotatable_bonds',
                                          'frac_csp3', 'defined_stereocentres')})
