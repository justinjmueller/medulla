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
        resp = input(f"[AUTH] Authenticate manually then press Enter, or 'n' to skip {experiment}: ")
        return resp.strip().lower() != 'n'

    cmd = [
        'htgettoken',
        '-a', 'htvaultprod.fnal.gov',
        '--vaulttokenttl=1d',
        '--vaulttokenminttl=12h',
        '-i', experiment,
    ]
    print(f"[AUTH] Running: {' '.join(cmd)}")
    # Do not capture output — htgettoken prints the OIDC URL the user must
    # visit, and waits for a response in the browser.  Suppressing it would
    # leave the process hanging with no visible prompt.
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"[AUTH] htgettoken failed for {experiment} (exit {result.returncode}).")
        return False

    print(f"[AUTH] Successfully authenticated for {experiment}.")
    return True
