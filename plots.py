"""Make figures from results.json (run geobell.py first)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import geobell as gb

d = json.load(open("results.json"))
main, sweep = d["main"], d["sweep"]
C = {"ideal": "#9aa5b1", "noisy": "#d1495b", "protected+noisy": "#2e86ab"}

# 1. fidelity per layout
layouts = list(dict.fromkeys(r["layout"] for r in main))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
w = 0.26
for i, case in enumerate(C):
    rs = [next(r for r in main if r["layout"] == L and r["case"] == case) for L in layouts]
    xs = [j + (i - 1) * w for j in range(len(layouts))]
    ax[0].bar(xs, [r["fidelity"] for r in rs], w, yerr=[r["fidelity_err"] for r in rs],
              color=C[case], label=case, capsize=3)
    if case != "ideal":
        ax[1].bar(xs, [r["acceptance"] for r in rs], w, color=C[case], label=case)
labels = [f"{L.split(' (')[0]}\n({next(r for r in main if r['layout']==L)['port_distance']} hops)"
          for L in layouts]
for a in ax:
    a.set_xticks(range(len(layouts)), labels, fontsize=9)
    a.spines[["top", "right"]].set_visible(False)
ax[0].set_ylim(0.5, 1.02); ax[0].set_ylabel("Bell fidelity F"); ax[0].legend(frameon=False, fontsize=8)
ax[0].axhline(0.5, color="k", lw=0.5)
ax[0].set_title("Fidelity between diameter ports (~50-qubit chips)")
ax[1].set_ylim(0, 1.05); ax[1].set_ylabel("accepted shot fraction")
ax[1].set_title("Cost of protection: postselection acceptance")
fig.tight_layout(); fig.savefig("fig_fidelity.png", dpi=150)

# 2. distance sweep
fig, ax = plt.subplots(figsize=(6.5, 4))
for case in ("noisy", "protected+noisy"):
    rs = [r for r in sweep if r["case"] == case]
    ax.errorbar([r["port_distance"] for r in rs], [r["fidelity"] for r in rs],
                yerr=[r["fidelity_err"] for r in rs], marker="o", color=C[case], label=case, capsize=2)
rs = [r for r in sweep if r["case"] == "protected+noisy"]
ax2 = ax.twinx()
ax2.plot([r["port_distance"] for r in rs], [r["acceptance"] for r in rs], "--", color="#2e86ab", alpha=0.5)
ax2.set_ylabel("acceptance (dashed)", color="#2e86ab"); ax2.set_ylim(0, 1.05)
for L in layouts:
    r = next(r for r in main if r["layout"] == L)
    ax.axvline(r["port_distance"], color="k", lw=0.6, ls=":")
    ax.text(r["port_distance"], 0.52, L.split(" (")[0].replace("Contemporary IBM", "heavy-hex").replace("2016-inspired", "ladder"),
            rotation=90, fontsize=7, ha="right", va="bottom")
ax.set_xlabel("graph distance between ports (hops)"); ax.set_ylabel("Bell fidelity F")
ax.set_ylim(0.5, 1.01); ax.legend(frameon=False, loc="upper right", fontsize=8)
ax.set_title("Fidelity vs. port distance (heavy-hex, same noise)")
ax.spines[["top"]].set_visible(False)
fig.tight_layout(); fig.savefig("fig_distance_sweep.png", dpi=150)

# 3. the three graphs with ports, ancilla and shortest path
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for a, (L, build) in zip(axs, gb.LAYOUTS.items()):
    g = build()
    p = gb.diameter_ports(g); anc = gb.geometry_aware_ancilla(g, *p)
    path = nx.shortest_path(g, *p)
    pos = (nx.kamada_kawai_layout(g) if "hyper" in L.lower() else
           {n: (n % 24, -(n // 24)) for n in g} if "ladder" in L else nx.kamada_kawai_layout(g))
    nx.draw_networkx_edges(g, pos, ax=a, edge_color="#ccc")
    nx.draw_networkx_edges(g, pos, edgelist=list(zip(path, path[1:])), ax=a, edge_color="#d1495b", width=2.5)
    col = ["#d1495b" if n in p else "#2e86ab" if n == anc else "#555" for n in g]
    size = [90 if n in p or n == anc else 18 for n in g]
    nx.draw_networkx_nodes(g, pos, ax=a, node_color=col, node_size=size)
    a.set_title(f"{L}\n{g.number_of_nodes()} qubits, ports {p} at {len(path)-1} hops", fontsize=10)
    a.axis("off")
fig.text(0.5, 0.02, "red = ports and routing path,  blue = geometry-aware ancilla (min total distance to both ports)",
         ha="center", fontsize=9)
fig.tight_layout(); fig.savefig("fig_layouts.png", dpi=150)
print("figures written")
