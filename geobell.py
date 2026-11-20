"""
GeoBell: geometry-aware Bell pairs across processor graphs.

WQCC Qiskit Fall Fest 2026, quanTA / USask prompt.

For three processor graphs (2016-style ladder, modern heavy-hex, hyperbolic {8,3}
lattice) we create a Bell pair between two fixed "port" qubits and compare:

  * ideal     : routed circuit, noiseless simulation
  * noisy     : routed circuit, user-defined NoiseModel
  * protected : routed circuit + XX and ZZ stabilizer checks measured by an ancilla
                placed at a geometry-aware position, followed by syndrome postselection

Fidelity follows the prompt: F = [1 + <XX> - <YY> + <ZZ>] / 4, where
<PP> = (agree - disagree) / shots, computed on accepted shots only when postselecting.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict

import networkx as nx
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
from qiskit.transpiler import CouplingMap
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError

SEED = 11
SHOTS = 8000
BASIS = ["cx", "h", "s", "sdg", "x", "y", "z", "reset", "measure"]  # Clifford basis -> stabilizer sim


# --------------------------------------------------------------------------- graphs
def ladder_2016(n_cols: int = 24) -> nx.Graph:
    """2 x N ladder, the layout style of the 16-qubit IBM QX5 (2x8), scaled to 48 qubits."""
    g = nx.Graph()
    for c in range(n_cols):
        g.add_edge(c, c + n_cols)                      # rungs
        if c < n_cols - 1:
            g.add_edge(c, c + 1)                       # top rail
            g.add_edge(c + n_cols, c + n_cols + 1)     # bottom rail
    return g


def heavy_hex(distance: int = 5) -> nx.Graph:
    """Heavy-hex lattice used by contemporary IBM processors (Eagle/Heron family)."""
    cm = CouplingMap.from_heavy_hex(distance, bidirectional=False)
    g = nx.Graph()
    g.add_edges_from(cm.get_edges())
    return g


def hyperbolic_83() -> nx.Graph:
    """A two-ring patch of the {8,3} hyperbolic tiling (octagons, 3 per vertex).

    Centre octagon c0..c7; each centre vertex gets one outward edge to o_i; each centre
    edge (c_i, c_{i+1}) is the inner edge of an outer octagon
        c_i - c_{i+1} - o_{i+1} - n1 - n2 - n3 - n4 - o_i - c_i.
    Three octagons meet at every interior vertex (3 x 135 deg = 405 deg > 360 deg),
    which is what makes the patch negatively curved. 48 qubits.
    """
    g = nx.Graph()
    c = list(range(8))
    o = list(range(8, 16))
    nxt = 16
    for i in range(8):
        g.add_edge(c[i], c[(i + 1) % 8])
        g.add_edge(c[i], o[i])
    for i in range(8):
        n = list(range(nxt, nxt + 4))
        nxt += 4
        path = [o[(i + 1) % 8], *n, o[i]]
        nx.add_path(g, path)
    return g


def diameter_ports(g: nx.Graph) -> tuple[int, int]:
    """Fixed ports at opposite 'ends' of the chip: the lowest-index diameter pair."""
    return pick_ports(g, nx.diameter(g))


def pick_ports(g: nx.Graph, target_distance: int) -> tuple[int, int]:
    """Fixed ports: lowest-index pair at the requested graph distance (deterministic)."""
    dist = dict(nx.all_pairs_shortest_path_length(g))
    for a in sorted(g.nodes):
        for b in sorted(g.nodes):
            if a < b and dist[a][b] == target_distance:
                return a, b
    raise ValueError("no pair at that distance")


def geometry_aware_ancilla(g: nx.Graph, a: int, b: int) -> int:
    """Ancilla = qubit minimising (dist to port A + dist to port B), tie-break balanced."""
    da = nx.single_source_shortest_path_length(g, a)
    db = nx.single_source_shortest_path_length(g, b)
    cands = [q for q in g.nodes if q not in (a, b)]
    return min(cands, key=lambda q: (da[q] + db[q], abs(da[q] - db[q]), q))


# --------------------------------------------------------------------------- circuits
def bell_circuit(basis: str, protected: bool) -> QuantumCircuit:
    """Logical circuit. q0 = port A, q1 = port B, q2 = ancilla (protected only)."""
    q = QuantumRegister(3 if protected else 2, "q")
    out = ClassicalRegister(2, "out")
    regs = [q, out]
    if protected:
        syn = ClassicalRegister(2, "syn")
        regs.append(syn)
    qc = QuantumCircuit(*regs)
    qc.h(0)
    qc.cx(0, 1)                     # |Phi+> = (|00> + |11>)/sqrt2, needs routing A <-> B

    if protected:
        # |Phi+> is the +1 eigenstate of both ZZ and XX (they commute), so measuring them
        # does not disturb the ideal state but flags any single X, Y or Z error on a port.
        # ZZ parity check -> catches bit flips
        qc.cx(0, 2)
        qc.cx(1, 2)
        qc.measure(2, syn[0])
        qc.reset(2)
        # XX parity check -> catches phase flips
        qc.h(2)
        qc.cx(2, 0)
        qc.cx(2, 1)
        qc.h(2)
        qc.measure(2, syn[1])

    if basis == "XX":
        qc.h([0, 1])
    elif basis == "YY":
        qc.sdg([0, 1])
        qc.h([0, 1])
    qc.measure(0, out[0])
    qc.measure(1, out[1])
    return qc


# --------------------------------------------------------------------------- noise
def make_noise(p1: float = 1e-3, p2: float = 1e-2, p_ro: float = 0.015) -> NoiseModel:
    """User-defined noise, identical on every graph so geometry is the only variable."""
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(p1, 1), ["h", "s", "sdg", "x", "y", "z"])
    nm.add_all_qubit_quantum_error(depolarizing_error(p2, 2), ["cx"])
    nm.add_all_qubit_readout_error(ReadoutError([[1 - p_ro, p_ro], [p_ro, 1 - p_ro]]))
    return nm


# --------------------------------------------------------------------------- metrics
def correlator(counts: dict, postselect: bool) -> tuple[float, int]:
    agree = disagree = 0
    for key, n in counts.items():
        parts = key.split()
        out = parts[-1]                       # 'out' register is listed last in keys
        if postselect and parts[0] != "00":  # 'syn' register first
            continue
        if out[0] == out[1]:
            agree += n
        else:
            disagree += n
    kept = agree + disagree
    return ((agree - disagree) / kept if kept else float("nan")), kept


@dataclass
class Row:
    layout: str
    case: str
    fidelity: float
    fidelity_err: float
    depth: int
    two_qubit_gates: int
    swaps: int
    physical_qubits: int
    ancillas: int
    acceptance: float
    port_distance: int
    ports: tuple
    ancilla_qubit: int | None


def compile_counts(qc: QuantumCircuit, cmap: CouplingMap, layout: list[int]):
    """Route once with SWAP kept as a gate (to count SWAPs), then lower to the basis."""
    routed = transpile(qc, coupling_map=cmap, initial_layout=layout,
                       basis_gates=BASIS + ["swap"], optimization_level=1,
                       seed_transpiler=SEED)
    swaps = routed.count_ops().get("swap", 0)
    lowered = transpile(routed, coupling_map=cmap, basis_gates=BASIS,
                        optimization_level=1, seed_transpiler=SEED)
    return lowered, swaps


def run_case(name: str, g: nx.Graph, ports: tuple[int, int], case: str,
             noise: NoiseModel | None) -> Row:
    protected = case == "protected+noisy"
    cmap = CouplingMap([list(e) for e in g.edges] + [[b, a] for a, b in g.edges])
    anc = geometry_aware_ancilla(g, *ports) if protected else None
    layout = [ports[0], ports[1]] + ([anc] if protected else [])

    sim = AerSimulator(method="stabilizer", noise_model=noise, seed_simulator=SEED)
    corr, var, kept_total, metrics = {}, 0.0, 0, None
    for basis in ("ZZ", "XX", "YY"):
        tqc, swaps = compile_counts(bell_circuit(basis, protected), cmap, layout)
        if basis == "ZZ":
            ops = tqc.count_ops()
            active = {tqc.find_bit(q).index for inst in tqc.data for q in inst.qubits}
            metrics = dict(depth=tqc.depth(), two_qubit_gates=ops.get("cx", 0),
                           swaps=swaps, physical_qubits=len(active))
        counts = sim.run(tqc, shots=SHOTS).result().get_counts()
        e, kept = correlator(counts, postselect=protected)
        corr[basis] = e
        var += (1 - e * e) / max(kept, 1)        # binomial variance of a +/-1 mean
        kept_total += kept

    fid = (1 + corr["XX"] - corr["YY"] + corr["ZZ"]) / 4
    return Row(layout=name, case=case, fidelity=fid, fidelity_err=math.sqrt(var) / 4,
               ancillas=1 if protected else 0,
               acceptance=kept_total / (3 * SHOTS),
               port_distance=nx.shortest_path_length(g, *ports), ports=ports,
               ancilla_qubit=anc, **metrics)


LAYOUTS = {
    "2016-inspired (2x24 ladder)": ladder_2016,
    "Contemporary IBM (heavy-hex)": heavy_hex,
    "Hyperbolic {8,3} patch": hyperbolic_83,
}


def main(noise_scale: float = 1.0) -> list[Row]:
    noise = make_noise(1e-3 * noise_scale, 1e-2 * noise_scale, 0.015 * noise_scale)
    rows = []
    for name, build in LAYOUTS.items():
        g = build()
        ports = diameter_ports(g)
        for case, nm in (("ideal", None), ("noisy", noise), ("protected+noisy", noise)):
            rows.append(run_case(name, g, ports, case, nm))
    return rows


def distance_sweep(max_d: int = 16) -> list[Row]:
    """Same noise, heavy-hex graph, ports moved further apart one hop at a time."""
    noise = make_noise()
    g = heavy_hex()
    rows = []
    for d in range(1, max_d + 1):
        ports = pick_ports(g, d)
        for case in ("noisy", "protected+noisy"):
            rows.append(run_case(f"heavy-hex d={d}", g, ports, case, noise))
    return rows


if __name__ == "__main__":
    rows = main()
    print(f"{'layout':30} {'case':16} {'F':>14} {'hops':>4} {'depth':>5} {'2q':>4} "
          f"{'swap':>4} {'phys':>4} {'acc':>5}")
    for r in rows:
        print(f"{r.layout:30} {r.case:16} {r.fidelity:.4f} ± {r.fidelity_err:.4f} "
              f"{r.port_distance:4d} {r.depth:5d} {r.two_qubit_gates:4d} {r.swaps:4d} "
              f"{r.physical_qubits:4d} {r.acceptance:5.2f}")
    sweep = distance_sweep()
    with open("results.json", "w") as f:
        json.dump({"main": [asdict(r) for r in rows],
                   "sweep": [asdict(r) for r in sweep]}, f, indent=2)
    print("wrote results.json")
