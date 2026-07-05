#!/usr/bin/env python3
"""
Generate a TOML systematics template from a CAF flat ROOT file.

uproot cannot handle memberwise serialisation of std::vector<SRWeightMapEntry>,
so PyROOT is used to read the weight parameter sets directly from the
globalTree.
"""
from argparse import ArgumentParser
from pathlib import Path

import ROOT

# Mapping from CAF weight-pset type integer to TOML type string.
_TYPE_NAME = {0: 'multisim', 3: 'multisigma'}

# ANSI helpers (no third-party dependency)
_INFO  = '\033[1m\033[94m[INFO]\033[0m'   # bold blue
_WARN  = '\033[1m\033[93m[WARN]\033[0m'   # bold yellow
_ERROR = '\033[1m\033[91m[ERROR]\033[0m'  # bold red

_TOML_HEADER = """\
[input]
path = 'output.root'
weights = 'data/*flat*.root'
use_additional_hash = true

[output]
path = 'output_sys.root'
histogram_destination = 'variations/'

[[sysvar]]
name = "reco_visible_energy"
bins = [25, 0, 3]
"""


def _str(s) -> str:
    """Decode a ROOT/C++ string to Python str, falling back to latin-1."""
    if isinstance(s, (bytes, bytearray)):
        return s.decode('latin-1')
    try:
        return str(s)
    except UnicodeDecodeError:
        return s.encode('latin-1', errors='replace').decode('latin-1')


def _snap_float(v: float, rtol: float = 1e-5, max_decimals: int = 10) -> float:
    """
    Snap a float to the shortest decimal representation within rtol of v.

    CAF ROOT branches store weights as 32-bit floats. When cppyy widens them
    to Python float (64-bit), values like 0.1f become 0.10000000149011612.
    This function detects such artefacts and returns the clean decimal (0.1)
    by trying round(v, d) for increasing d and accepting the first result
    whose relative error falls below rtol.

    Parameters
    ----------
    v : float
        Value to snap.
    rtol : float
        Relative tolerance for accepting a rounded value (default: 1e-5).
    max_decimals : int
        Maximum number of decimal places to try (default: 10).

    Returns
    -------
    float
        Snapped value if a clean representation was found, otherwise v.
    """
    if v == 0.0:
        return v
    for d in range(0, max_decimals + 1):
        r = round(v, d)
        if r != 0.0 and abs(v - r) / abs(r) < rtol:
            return r
    return v


def read_weight_psets(input_file: str) -> list[dict]:
    """
    Open a CAF flat ROOT file and extract the weight parameter sets stored in
    the global tree as plain Python dicts.

    All C++ data is copied into native Python types before the file is closed,
    avoiding dangling-pointer issues that arise when cppyy wrappers outlive
    the owning TFile.

    Parameters
    ----------
    input_file : str
        Path (or XRootD URL) to the CAF flat ROOT file.

    Returns
    -------
    psets : list[dict]
        List of dicts, one per caf::SRWeightPSet, with keys:

        ``name`` (str)
            Weight parameter set name.
        ``type`` (int)
            Integer type code (0 = multisim, 3 = multisigma).
        ``nsigma`` (list[float], multisigma only)
            Sigma points collected across all map entries.
    """
    rf = ROOT.TFile.Open(input_file)
    if not rf or rf.IsZombie():
        raise FileNotFoundError(f"Cannot open ROOT file: {input_file}")
    tree = rf.Get('globalTree')
    if not tree:
        raise KeyError("'globalTree' not found in ROOT file.")
    tree.GetEntry(0)
    sr_global = getattr(tree, 'global')

    psets = []
    for pset in sr_global.wgts:
        entry: dict = {
            'name': _str(pset.name),
            'type': int(pset.type),
        }
        if int(pset.type) == 3:  # multisigma — collect nsigma points now
            nsigma: list[float] = []
            for map_entry in pset.map:
                nsigma.extend(_snap_float(float(v)) for v in map_entry.vals)
            entry['nsigma'] = nsigma
        psets.append(entry)

    rf.Close()
    return psets


def build_toml(psets: list[dict]) -> tuple[str, int]:
    """
    Build a TOML systematics template string from a list of parameter-set
    dicts (as returned by :func:`read_weight_psets`).

    Parameters
    ----------
    psets : list[dict]
        List of parameter-set dicts with at minimum ``name`` and ``type``
        keys, and an ``nsigma`` key for multisigma entries.

    Returns
    -------
    toml_str : str
        The complete TOML template as a string, ready to be written to disk.
    skipped : int
        Number of parameter sets skipped due to an unrecognised type.
    """
    blocks = [_TOML_HEADER]
    skipped = 0

    for i, pset in enumerate(psets):
        type_str = _TYPE_NAME.get(pset['type'])
        if type_str is None:
            print(
                f"{_WARN} -- Unknown type {pset['type']} for pset {i} "
                f"({pset['name']}), skipping."
            )
            skipped += 1
            continue

        block = [
            '[[sys]]',
            f'name = "{pset["name"]}"',
            f'type = "{type_str}"',
            f'index = {i}',
        ]

        if type_str == 'multisigma':
            block.append(f'nsigma = {pset["nsigma"]}')

        blocks.append('\n'.join(block))

    return '\n\n'.join(blocks) + '\n', skipped


def main(
    input_file: str,
    output_file: str,
) -> None:
    """
    Generate a TOML systematics template by reading weight parameter sets
    from a CAF flat ROOT file.

    Parameters
    ----------
    input_file : str
        Path (or XRootD URL) to the input CAF flat ROOT file.
    output_file : str
        Path to write the generated TOML template to.

    Returns
    -------
    None.
    """
    print(f"{_INFO} -- Reading weight parameter sets from: {input_file}")
    psets = read_weight_psets(input_file)

    toml_str, skipped = build_toml(psets)

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(toml_str)

    n_written = len(psets) - skipped
    suffix = f" ({skipped} skipped)." if skipped else "."
    print(f"{_INFO} -- Wrote {n_written} systematics to {out_path}{suffix}")


if __name__ == '__main__':
    p = ArgumentParser(
        description='Generate a TOML systematics template from a CAF flat ROOT file.'
    )
    p.add_argument(
        '--input', '-i', type=str, required=True,
        help='Path (or XRootD URL) to the input CAF flat ROOT file.',
    )
    p.add_argument(
        '--output', '-o', type=str, default='sys_template_generated.toml',
        help='Path for the generated TOML template (default: sys_template_generated.toml).',
    )
    args = p.parse_args()
    main(input_file=args.input, output_file=args.output)
