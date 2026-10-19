"""Tests for scripts/guard_bash.py (the PreToolUse Bash guard)."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import guard_bash  # noqa: E402

P = '/pnfs/icarus/scratch/users/micarrig/nueCCInclusiveSpine31'


@pytest.mark.parametrize('cmd', [
    'pkill -f plotEfficiencyCLI.py',
    'killall python3',
    'source ~/scripts/setup_icarus_env.sh',
    '. ~/scripts/setup_icarus_env.sh >/dev/null 2>&1; python3 x.py',
    'bash ~/scripts/setup_icarus_env.sh',
    'jobsub_submit -G icarus file:///x.sh',
    'jobsub_rm --jobid 1@x',
    f'rm -rf {P}/output',
    f'mv {P}/output /tmp/x',
    'rm output_jobid0001.root /pnfs/icarus/scratch/users/micarrig/foo/',
    f'find {P}/output -name "*.root" -delete',
    f'find {P} -size -1k -exec rm {{}} \\;',
])
def test_denied(cmd):
    assert guard_bash.decide(cmd)[0] == 'deny'


@pytest.mark.parametrize('cmd', [
    f'python3 batch/medulla.py -p {P} -e icarus --launch-jobs --tag abc',
    f'python3 batch/medulla.py -p {P} -e icarus -Z -y --variation-phase2 -V x.toml',
    f'python3 batch/medulla.py -p {P} -e icarus -c -t sel.toml -b 30',
    'python3 analysis/plotEfficiencyCLI.py --toml a.toml --input f.root --output-dir d',
    'python3 analysis/plotEfficiency.py',
    'python3 analysis/plot_all.py --only selection',
    'python3 analysis/makeSelectionPlotsCLI.py --toml a.toml --input f.root --output-dir d',
    'python3 analysis/makeNMinus1PlotsCLI.py --toml a.toml --input f.root --output-dir d',
    'python3 batch/merge_list.py /x --tags detector_variation',
    '~/scripts/mergeFiles.sh -d -l list.txt out.root',
    'scripts/in_sl7.sh "hadd -f out.root a.root b.root"',
    'python3 systematics/add_numi_flux_sys.py merged.root events/NuMIFull selected',
    'build/systematics/to_gundam cfg.toml merged.root',
    'export BEARER_TOKEN_FILE=/tmp/bt_u1 && htgettoken -a htvaultprod.fnal.gov -i icarus',
    'rm -rf /tmp/mydir',
    'mv old.toml new.toml',
    'find /tmp/x -name "*.tmp" -delete',
    'python3 -c "import shutil; shutil.rmtree(\'/tmp/x\')"',
    'python3 -c "import os; os.remove(\'a\')"',
])
def test_ask(cmd):
    assert guard_bash.decide(cmd)[0] == 'ask'


@pytest.mark.parametrize('cmd', [
    f'python3 batch/medulla.py -p {P} -e icarus --launch-jobs --dry-run',
    f'python3 batch/medulla.py -p {P} -e icarus -y --variation-phase2 -V x.toml --dry-run --tag abc',
    'python3 batch/medulla.py --help',
    'python3 batch/medulla.py -h',
    'python3 analysis/plotEfficiencyCLI.py --toml a.toml --input f.root --output-dir d --dry-run',
    'python3 analysis/plot_all.py --dry-run',
    'python3 analysis/makeSelectionPlotsCLI.py --toml a.toml --input f.root --output-dir d --dry-run',
    'python3 batch/merge_list.py /x --tags detector_variation --dry-run',
])
def test_dry_run_and_help_pass(cmd):
    assert guard_bash.decide(cmd) is None


@pytest.mark.parametrize('cmd', [
    'ls -la /pnfs/icarus/scratch/users/micarrig',
    'python3 -m pytest batch -q',
    'git status --short',
    'git ls-remote origin feature/carriganm95_working',
    'jobsub_q --user micarrig',
    'pgrep -fa "[p]lotEfficiencyCLI"',
    'free -g',
    'grep -n medulla CLAUDE.md',                 # mentions the name, not the script path
    'cat systematics/src/to_gundam.cc | head',   # a source file, not the binary
    'scripts/in_sl7.sh "which root; root-config --version"',
])
def test_ordinary_commands_pass(cmd):
    assert guard_bash.decide(cmd) is None


@pytest.mark.parametrize('cmd', [
    'cat batch/medulla.py | head -5',   # conservative: a path mention asks; harmless
])
def test_path_mentions_ask_conservatively(cmd):
    assert guard_bash.decide(cmd)[0] == 'ask'


def run_hook(payload):
    r = subprocess.run([sys.executable, str(SCRIPTS / 'guard_bash.py')], input=payload,
                       capture_output=True, text=True)
    assert r.returncode == 0
    return r.stdout


def test_hook_output_format_deny():
    out = json.loads(run_hook(json.dumps({'tool_name': 'Bash', 'tool_input': {'command': 'pkill -f x'}})))
    h = out['hookSpecificOutput']
    assert h['hookEventName'] == 'PreToolUse'
    assert h['permissionDecision'] == 'deny'
    assert h['permissionDecisionReason'].startswith('[guard_bash]')


def test_hook_never_allows():
    for cmd in ['ls', 'pkill x', 'python3 batch/medulla.py -c', 'python3 batch/medulla.py --dry-run']:
        out = run_hook(json.dumps({'tool_name': 'Bash', 'tool_input': {'command': cmd}}))
        if out:
            assert json.loads(out)['hookSpecificOutput']['permissionDecision'] in ('ask', 'deny')


@pytest.mark.parametrize('payload', ['', 'not json', '{}', '{"tool_name": "Read"}',
                                     '{"tool_name": "Bash", "tool_input": {"command": 5}}'])
def test_hook_fails_open(payload):
    assert run_hook(payload) == ''


def test_kill_asks():
    assert guard_bash.decide('kill -TERM 51590 51587')[0] == 'ask'
    assert guard_bash.decide('kill 123')[0] == 'ask'
    assert guard_bash.decide('echo skill') is None


def _decision(cmd, mode):
    payload = {'tool_name': 'Bash', 'tool_input': {'command': cmd}}
    if mode is not None:
        payload['permission_mode'] = mode
    out = run_hook(json.dumps(payload))
    return json.loads(out)['hookSpecificOutput']['permissionDecision'] if out else None


@pytest.mark.parametrize('mode', ['auto', 'bypassPermissions', 'dontAsk'])
def test_ask_becomes_deny_when_nobody_is_prompted(mode):
    for cmd in ['kill 123', 'python3 batch/medulla.py -p x -e icarus --launch-jobs',
                'python3 batch/merge_list.py P --output m.txt',
                # plotting chained with something strict is still blocked
                'python3 analysis/plot_all.py --only nminus1; kill 123']:
        assert _decision(cmd, mode) == 'deny'
    # Dry-runs and ordinary commands are untouched.
    assert _decision('python3 analysis/plot_all.py --dry-run', mode) is None
    assert _decision('ls', mode) is None


def test_plotting_left_to_classifier_in_auto_mode_only():
    cmd = 'python3 analysis/plot_all.py --only nminus1'
    assert _decision(cmd, 'auto') == 'ask'
    assert _decision(cmd, 'bypassPermissions') == 'deny'
    assert _decision(cmd, 'dontAsk') == 'deny'


@pytest.mark.parametrize('mode', [None, 'default', 'acceptEdits', 'plan'])
def test_ask_stays_ask_when_owner_is_prompted(mode):
    assert _decision('python3 analysis/plot_all.py --only nminus1', mode) == 'ask'
    assert _decision('kill 123', mode) == 'ask'
