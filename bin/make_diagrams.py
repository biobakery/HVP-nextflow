#!/usr/bin/env python3
"""Generate the HVP architecture diagram and component table from the code.

Nothing here is hand-drawn. Two sources are read:

  conf/components.json   the registry: stages, the data each one reads and
                         writes, and the components in each (with status)
  workflows/hvp.nf       the orchestrator, to confirm every registered
                         component is actually included and wired

and three things are written:

  docs/diagrams/architecture.mmd       the Mermaid source
  README.md, docs/architecture.md      the diagram and the component table,
                                       between <!-- generated:... --> markers

The diagram draws one arrow per hand-off of data between stages, labelled with
what is handed over, instead of one per process. Components inside a "chain"
stage get an arrow between them because they run in that order; components in
a "parallel" stage get none because they don't depend on each other.

    bin/make_diagrams.py           regenerate
    bin/make_diagrams.py --check   fail if anything is stale or out of sync (CI)
    bin/make_diagrams.py --render  also write architecture.svg/.png (needs mermaid-cli, `mmdc`)

Needs only Python 3: no Nextflow, no tools. --render additionally needs mmdc on
PATH (npm i -g @mermaid-js/mermaid-cli); set PUPPETEER_CONFIG to a puppeteer
JSON config to point it at an installed Chrome.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REGISTRY = os.path.join(REPO, "conf", "components.json")
ORCHESTRATOR = os.path.join(REPO, "workflows", "hvp.nf")
COMPONENTS_DIR = os.path.join(REPO, "subworkflows", "local", "components")
MMD = os.path.join(REPO, "docs", "diagrams", "architecture.mmd")
PAGES = [os.path.join(REPO, "README.md"), os.path.join(REPO, "docs", "architecture.md")]


# ── Consistency: registry ⇄ component files ⇄ orchestrator ───────────────────

def consistency_errors(reg):
    """Everything that would make the diagram lie about the code."""
    errors = []
    orch = open(ORCHESTRATOR).read()
    stage_ids = [s["id"] for s in reg["stages"]]
    data = reg["data"]

    for s in reg["stages"]:
        for key in ("consumes", "produces"):
            if s[key] not in data:
                errors.append(f"stage {s['id']}: {key} '{s[key]}' is not in data")
        if s["default"] and not any(c["name"] == s["default"] and c["stage"] == s["id"]
                                    for c in reg["components"]):
            errors.append(f"stage {s['id']}: default '{s['default']}' is not one of its components")

    for c in reg["components"]:
        name, wf = c["name"], c["name"].upper()
        if c["stage"] not in stage_ids:
            errors.append(f"{name}: unknown stage '{c['stage']}'")
        if c["status"] not in ("ready", "skeleton"):
            errors.append(f"{name}: status must be ready or skeleton, not '{c['status']}'")
        path = os.path.join(COMPONENTS_DIR, name, "main.nf")
        if not os.path.exists(path):
            errors.append(f"{name}: missing {os.path.relpath(path, REPO)}")
            continue
        text = open(path).read()
        if not re.search(rf"^workflow\s+{wf}\s*\{{", text, re.M):
            errors.append(f"{name}: {os.path.relpath(path, REPO)} does not define workflow {wf}")
        header = re.search(r"//\s*Component:\s*(\w+)\s*·\s*stage:\s*(\w+)\s*·\s*status:\s*(\w+)", text)
        if not header:
            errors.append(f"{name}: no '// Component: … · stage: … · status: …' header line")
        elif header.groups() != (name, c["stage"], c["status"]):
            errors.append(f"{name}: header says {header.groups()}, registry says "
                          f"({name}, {c['stage']}, {c['status']})")
        if not re.search(rf"include\s*\{{\s*{wf}\s*\}}\s*from\s*'[^']*/components/{name}/main\.nf'", orch):
            errors.append(f"{name}: not included in workflows/hvp.nf")
        if f"on('{name}')" not in orch:
            errors.append(f"{name}: no on('{name}') block in workflows/hvp.nf")

    registered = {c["name"] for c in reg["components"]}
    for d in sorted(os.listdir(COMPONENTS_DIR)):
        if os.path.isdir(os.path.join(COMPONENTS_DIR, d)) and d not in registered:
            errors.append(f"subworkflows/local/components/{d} is not in conf/components.json")
    return errors


# ── Diagram ──────────────────────────────────────────────────────────────────

def mermaid(reg):
    stages = reg["stages"]
    data = reg["data"]
    by_stage = {s["id"]: [c for c in reg["components"] if c["stage"] == s["id"]] for s in stages}
    producer = {s["produces"]: s["id"] for s in stages}
    consumed = {s["consumes"] for s in stages}

    out = [
        '%%{init: {"flowchart": {"curve": "basis", "nodeSpacing": 24, "rankSpacing": 44, "padding": 10}}}%%',
        "flowchart TB",
        '    entry(["<b>main.nf</b> → orchestrator<br/><i>--components</i> · <i>--plugin</i>"])',
    ]

    for i, s in enumerate(stages, 1):
        how = "in order" if s["mode"] == "chain" else "side by side"
        out.append(f'    subgraph st_{s["id"]}["<b>{i} · {s["title"]}</b> &nbsp;<i>{how}</i>"]')
        out.append("        direction LR")
        for c in by_stage[s["id"]]:
            tools = ", ".join(c["tools"])
            out.append(f'        {c["name"]}["<b>{c["name"]}</b><br/><i>{tools}</i>"]')
        if s["mode"] == "chain":
            names = [c["name"] for c in by_stage[s["id"]]]
            for a, b in zip(names, names[1:]):
                out.append(f"        {a} --> {b}")
        out.append("    end")

    # One arrow per hand-off between stages, labelled with the data.
    for s in stages:
        src = "entry" if s["consumes"] == "raw_reads" else "st_" + producer[s["consumes"]]
        out.append(f'    {src} -- "{data[s["consumes"]]}" --> st_{s["id"]}')

    # Data nothing downstream reads is a pipeline output.
    for s in stages:
        if s["produces"] not in consumed:
            out.append(f'    st_{s["id"]} --> out_{s["produces"]}[("{data[s["produces"]]}")]')

    out += [
        '    subgraph legend[" "]',
        "        direction LR",
        '        l_ready["ready"]',
        '        l_skel["placeholder"]',
        "    end",
        "",
        "    classDef ready fill:#dbeafe,stroke:#2563eb,stroke-width:1.5px,color:#0f172a",
        "    classDef skeleton fill:#ffffff,stroke:#94a3b8,stroke-width:1px,stroke-dasharray:4 3,color:#64748b",
        "    classDef entry fill:#e2e8f0,stroke:#475569,stroke-width:1.5px,color:#0f172a",
        "    classDef output fill:#dcfce7,stroke:#16a34a,stroke-width:1.5px,color:#14532d",
    ]
    ready = [c["name"] for c in reg["components"] if c["status"] == "ready"] + ["l_ready"]
    skel = [c["name"] for c in reg["components"] if c["status"] == "skeleton"] + ["l_skel"]
    outs = [f"out_{s['produces']}" for s in stages if s["produces"] not in consumed]
    out.append(f"    class {','.join(ready)} ready")
    out.append(f"    class {','.join(skel)} skeleton")
    out.append("    class entry entry")
    out.append(f"    class {','.join(outs)} output")
    for s in stages:
        out.append(f"    style st_{s['id']} fill:#f8fafc,stroke:#cbd5e1,stroke-width:1px,color:#334155")
    out.append("    style legend fill:none,stroke:none")
    out.append("    linkStyle default stroke:#64748b,stroke-width:1.5px")
    return "\n".join(out) + "\n"


def component_table(reg):
    titles = {s["id"]: f"{i} · {s['title']}" for i, s in enumerate(reg["stages"], 1)}
    rows = ["| Stage | Component | Status | Tools | What it does |",
            "|---|---|---|---|---|"]
    for c in reg["components"]:
        status = "ready" if c["status"] == "ready" else "placeholder"
        rows.append(f"| {titles[c['stage']]} | `{c['name']}` | {status} | "
                    f"{', '.join(c['tools'])} | {c['summary']} |")
    return "\n".join(rows) + "\n"


# ── Writing ──────────────────────────────────────────────────────────────────

def splice(text, marker, body):
    """Replace what sits between <!-- generated:marker --> and <!-- /generated:marker -->."""
    start, end = f"<!-- generated:{marker} -->", f"<!-- /generated:{marker} -->"
    pattern = re.compile(re.escape(start) + r".*?" + re.escape(end), re.S)
    if not pattern.search(text):
        return text
    note = f"<!-- written by bin/make_diagrams.py from conf/components.json; do not edit by hand -->\n"
    return pattern.sub(lambda _: f"{start}\n{note}{body}{end}", text)


def render(reg):
    """Return {path: content} for everything this script owns."""
    diagram = mermaid(reg)
    table = component_table(reg)
    files = {MMD: diagram}
    for page in PAGES:
        if not os.path.exists(page):
            continue
        text = open(page).read()
        text = splice(text, "architecture", f"```mermaid\n{diagram}```\n")
        text = splice(text, "components", table)
        files[page] = text
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="fail if anything is stale or out of sync")
    ap.add_argument("--render", action="store_true", help="also render docs/diagrams/architecture.{svg,png}")
    args = ap.parse_args()

    reg = json.load(open(REGISTRY))
    errors = consistency_errors(reg)
    if errors:
        print("conf/components.json and the code disagree:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1

    stale = []
    for path, content in render(reg).items():
        current = open(path).read() if os.path.exists(path) else None
        if current == content:
            continue
        if args.check:
            stale.append(os.path.relpath(path, REPO))
        else:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w").write(content)
            print(f"wrote {os.path.relpath(path, REPO)}")

    if stale:
        print("Out of date; run bin/make_diagrams.py and commit:\n  " + "\n  ".join(stale), file=sys.stderr)
        return 1
    if args.check:
        print("diagram, component table and registry are in sync")
    if args.render:
        return render_images()
    return 0


def render_images():
    mmdc = shutil.which("mmdc")
    if not mmdc:
        print("--render needs mermaid-cli: npm i -g @mermaid-js/mermaid-cli", file=sys.stderr)
        return 1
    extra = ["-p", os.environ["PUPPETEER_CONFIG"]] if os.environ.get("PUPPETEER_CONFIG") else []
    for ext, opts in (("svg", []), ("png", ["-s", "2"])):
        out = MMD[:-len(".mmd")] + "." + ext
        subprocess.run([mmdc, "-i", MMD, "-o", out, "-b", "white", *opts, *extra], check=True)
        print(f"wrote {os.path.relpath(out, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
