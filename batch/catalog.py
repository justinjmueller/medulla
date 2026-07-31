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
    default_catalog_path = Path(catalog_path) if catalog_path is not None else DEFAULT_CATALOG

    # Cache loaded catalogs by resolved path, since multiple directives may
    # reference the same file.
    catalog_cache = {}

    def _load_catalog(path):
        path = Path(path)
        if path not in catalog_cache:
            catalog = toml.load(path)
            catalog_cache[path] = {
                s['key']: s for s in catalog.get('sample', []) if 'key' in s
            }
        return catalog_cache[path]

    for inc in includes:
        requested_keys = set(inc.get('keys', []))

        # Empty keys list → resolve nothing for this directive
        if not requested_keys:
            continue

        # Each directive may name its own catalog file (relative to the
        # shared catalog directory), falling back to the caller-supplied
        # default when omitted.
        inc_file = inc.get('file')
        if inc_file:
            this_catalog_path = Path(inc_file)
            if not this_catalog_path.is_absolute():
                this_catalog_path = CATALOG_DIR / inc_file
        else:
            this_catalog_path = default_catalog_path
        catalog_by_key = _load_catalog(this_catalog_path)

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
