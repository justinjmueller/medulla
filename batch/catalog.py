"""Sample catalog resolution for medulla campaign and batch tools."""
import toml
from pathlib import Path

CATALOG_DIR = Path(__file__).resolve().parent.parent / 'selection' / 'toml'

def resolve_samples(cfg, toml_dir=None, enable_keys=None):
    """
    Resolve [[include_samples]] directives in a parsed TOML config.

    Parameters
    ----------
    cfg : dict
        Parsed TOML configuration dictionary.
    toml_dir : Path
        Directory of the source TOML (for resolving relative paths).
    enable_keys : list[str]
        Sample keys to enable (disable=False). All others disabled.
        If None, uses the 'enable' list from the directive itself.

    Returns
    -------
    dict
        Config with [[include_samples]] replaced by [[sample]].
    """
    includes = cfg.pop('include_samples', [])
    if not includes:
        return cfg

    resolved = cfg.get('sample', [])
    # The catalog path is relative to the parent of the analysis directory
    # (i.e., selection/toml/), not the analysis directory itself.
    base = (toml_dir.parent if toml_dir else CATALOG_DIR)

    for inc in includes:
        catalog = toml.load(base / inc['file'])
        requested_keys = set(inc.get('keys', []))

        # Priority: function argument > directive 'enable' field > none
        if enable_keys is not None:
            active_keys = set(enable_keys)
        else:
            active_keys = set(inc.get('enable', []))

        for sample in catalog.get('sample', []):
            key = sample.get('key')
            if key is None:
                continue
            if requested_keys and key not in requested_keys:
                continue

            entry = {
                'name': sample['name'],
                'path': sample['path'],
                'ismc': sample['ismc'],
                'disable': key not in active_keys if active_keys else True,
            }
            resolved.append(entry)

    cfg['sample'] = resolved
    return cfg
