"""Authentication helper for medulla campaign management."""
import subprocess
import shutil


def authenticate(experiment):
    """
    Authenticate for the given experiment using htgettoken.
    Falls back to interactive prompt if htgettoken is unavailable or fails.

    Parameters
    ----------
    experiment : str
        The experiment to authenticate for ('sbnd' or 'icarus').

    Returns
    -------
    bool
        True if authentication succeeded or user confirmed manual auth.
        False if user chose to skip.
    """
    htgettoken = shutil.which('htgettoken')
    if htgettoken is None:
        print(f"[AUTH] htgettoken not found. Please authenticate manually for {experiment}.")
        input(f"[AUTH] Press Enter after authenticating for {experiment}...")
        return True

    cmd = [
        'htgettoken',
        '-a', 'htvaultprod.fnal.gov',
        '--vaulttokenttl=1d',
        '--vaulttokenminttl=12h',
        '-i', experiment,
    ]
    print(f"[AUTH] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[AUTH] htgettoken failed for {experiment}: {result.stderr.strip()}")
        resp = input(f"[AUTH] Authenticate manually and press Enter, or 'skip' to skip {experiment}: ")
        return resp.strip().lower() != 'skip'

    print(f"[AUTH] Successfully authenticated for {experiment}.")
    return True
