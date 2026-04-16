"""Sample catalog resolution for medulla campaign and batch tools."""
import toml
from pathlib import Path

CATALOG_DIR = Path(__file__).resolve().parent.parent / 'selection' / 'toml'
DEFAULT_CATALOG = CATALOG_DIR / 'common' / 'samples.toml'

def resolve_samples(cfg, catalog_path=None, enable_keys=None):
    """
    Resolve [[include_samples]] directives in a parsed TOML config.

    Parameters
    ----------
    cfg : dict
        Parsed TOML configuration dictionary.
    catalog_path : str | Path
        Path to the shared sample catalog TOML file.
    enable_keys : list[str]
        Sample keys to enable (disable=False). All others disabled.
        If None, uses the 'enable' list from the directive itself.

    Returns
    -------
    dict
        Config with [[include_samples]] replaced by [[sample]].

    Raises
    ------
    KeyError
        If a requested key is not found in the catalog.
    FileNotFoundError
        If the catalog file does not exist.
    """
    includes = cfg.pop('include_samples', [])
    if not includes:
        return cfg

    resolved = cfg.get('sample', [])

    if catalog_path is None:
        catalog_path = DEFAULT_CATALOG
    catalog = toml.load(Path(catalog_path))

    catalog_by_key = {
        s['key']: s for s in catalog.get('sample', []) if 'key' in s
    }

    for inc in includes:
        requested_keys = set(inc.get('keys', []))

        # Empty keys list → resolve nothing for this directive
        if not requested_keys:
            continue

        # Priority: function argument > directive 'enable' field > none
        if enable_keys is not None:
            active_keys = set(enable_keys)
        else:
            active_keys = set(inc.get('enable', []))

        for key in requested_keys:
            if key not in catalog_by_key:
                raise KeyError(key)
            sample = catalog_by_key[key]
            entry = {
                'name': sample['name'],
                'path': sample['path'],
                'ismc': sample['ismc'],
                'disable': (not active_keys) or (key not in active_keys),
            }
            resolved.append(entry)

    cfg['sample'] = resolved
    return cfg
