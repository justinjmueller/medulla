"""Tests for scripts/guard_bash.py (the PreToolUse Bash guard)."""
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import guard_bash  # noqa: E402

P = '/pnfs/icarus/scratch/users/micarrig/nueCCInclusiveSpine31'
T = 'systematics/toml/NuMI_nue_gundam.toml'


@pytest.fixture
def big_node(monkeypatch):
    monkeypatch.setattr(guard_bash, 'available_gb', lambda: 344.0)


@pytest.fixture
def small_node(monkeypatch):
    monkeypatch.setattr(guard_bash, 'available_gb', lambda: 9.0)


def decide(cmd, **kw):
    kw.setdefault('allow', [])
    return guard_bash.decide(cmd, **kw)


# --- deny ------------------------------------------------------------------

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
    'echo x >> scripts/guard_allow.txt',
    'tee -a scripts/guard_allow.txt < x',
    "sed -i 's/a/b/' scripts/guard_allow.txt",
    'echo ok; pkill -f x',
])
def test_denied(cmd):
    assert decide(cmd)[0] == 'deny'


# --- ask: only real submissions, --yes, and the remaining hard cases --------

@pytest.mark.parametrize('cmd', [
    f'python3 batch/medulla.py -p {P} -e icarus -Z -y --launch-jobs --tag abc --memory 10000 --disk 10',
    f'python3 batch/medulla.py -p {P} -e icarus --launch-jobs 5 --tag abc',
    f'python3 batch/medulla.py -p {P} -e icarus -Z -y --variation-phase2 -V {T} --tag abc',
    f'python3 batch/medulla.py -p {P} -e icarus --variation-phase1 --test-job -V {T} --tag abc',  # grid phase 1
    f'python3 batch/medulla.py -p {P} -e icarus --variation-phase1 -V {T} --tag abc',   # grid phase 1
    f'python3 batch/medulla.py -p {P} -e icarus --launch-jobs -r --tag abc',
    f'python3 batch/medulla.py -p {P} -e icarus -Zr --tag abc',
    f'BEARER_TOKEN_FILE=/tmp/bt_u1 python3 batch/medulla.py -p {P} -e icarus --launch-jobs --tag abc',
])
def test_grid_submission_asks(cmd):
    d, reason = decide(cmd)
    assert d == 'ask' and 'grid jobs' in reason


@pytest.mark.parametrize('cmd', [
    f'python3 batch/medulla.py -p {P} -e icarus -Z -y < /dev/null',
    f'python3 batch/medulla.py -p {P} -e icarus --yes',
    f'python3 batch/medulla.py -p {P} -e icarus -Zy',
])
def test_medulla_yes_asks(cmd):
    d, reason = decide(cmd)
    assert d == 'ask' and '--yes' in reason


@pytest.mark.parametrize('cmd', [
    'python3 systematics/add_numi_flux_sys.py merged.root events/NuMIFull selected',
    "scripts/in_sl7.sh 'python3 systematics/add_numi_flux_sys.py m.root events/g4 selected'",
    'export BEARER_TOKEN_FILE=/tmp/bt_u1 && htgettoken -a htvaultprod.fnal.gov -i icarus',
    'rm -rf /tmp/mydir',
    'mv old.toml new.toml',
    'find /tmp/x -name "*.tmp" -delete',
    'python3 -c "import shutil; shutil.rmtree(\'/tmp/x\')"',
    'python3 -c "import os; os.remove(\'a\')"',
    'kill -TERM 51590 51587',
    'kill 123',
])
def test_still_asks(cmd):
    assert decide(cmd)[0] == 'ask'


# --- automatic: dry-runs, status checks, local stages, merges, plots --------

@pytest.mark.parametrize('cmd', [
    # every grid stage as a dry run, including the -y the real command carries
    f'python3 batch/medulla.py -p {P} -e icarus -Z -y --launch-jobs --tag abc --memory 10000 --disk 10 --dry-run',
    f'python3 batch/medulla.py -p {P} -e icarus -y --variation-phase2 -V {T} --dry-run --tag abc',
    f'python3 batch/medulla.py -p {P} -e icarus --test-job --variation-phase2 -V {T} --dry-run',
    f'export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u) && python3 batch/medulla.py -p {P} -e icarus --launch-jobs --dry-run 2>&1 < /dev/null | tee -a {P}/logs/stage2.log',
    'python3 batch/medulla.py --help',
    'python3 batch/medulla.py -h',
    # one test job (selection or phase 2), with the -y it needs for its prompt
    f'python3 batch/medulla.py -p {P} -e icarus -Z -y --test-job --tag abc --memory 10000 --disk 10',
    f'python3 batch/medulla.py -p {P} -e icarus -y --test-job --variation-phase2 -V {T} --tag abc',
    f'python3 batch/medulla.py -p {P} -e icarus -ZyT --tag abc',
    f'python3 batch/medulla.py -p {P} -e icarus -y --test-job -r --tag abc',
    # read-only status check, project creation, interactive phase 1
    f'python3 batch/medulla.py -p {P} -e icarus -Z < /dev/null',
    f'python3 batch/medulla.py -p {P} -e icarus -Z 2>&1 < /dev/null | tee -a {P}/logs/status.log',
    f'python3 batch/medulla.py -p {P} -e icarus -c -t sel.toml -s {T} -b 30',
    f"scripts/in_sl7.sh 'python3 batch/medulla.py -p {P} -e icarus -y --variation-phase1 --variation-interactive -V {T} -I {P}/merged_stage0.root'",
    # merges into new files, merge lists, stage 5 with --output, stage 6
    f"scripts/in_sl7.sh 'python3 batch/merge_list.py {P} --tags detector_variation --no-varsys --output {P}/merge_list_stage0.txt && ~/scripts/mergeFiles.sh -d -l {P}/merge_list_stage0.txt /tmp/does_not_exist/merged_stage0.root'",
    'python3 batch/merge_list.py /x --tags detector_variation',
    '~/scripts/mergeFiles.sh -d -l list.txt /tmp/does_not_exist/out.root',
    'scripts/in_sl7.sh "hadd /tmp/does_not_exist/out.root a.root b.root"',
    f"scripts/in_sl7.sh 'python3 systematics/add_numi_flux_sys.py {P}/merged_final.root --output {P}/merged_final_flux.root events/NuMIFull all_signal selected true_signal nc_sideband'",
    'python3 systematics/add_numi_flux_sys.py m.root --output=o.root events/g4 selected',
    f"scripts/in_sl7.sh 'mkdir -p {P}/gundam_tmp && cd {P}/gundam_tmp && /x/build/systematics/to_gundam /x/{T} {P}/merged_final_flux.root'",
    'build/systematics/to_gundam cfg.toml merged.root',
    # plot_all.py has its own memory gate: always automatic
    'python3 analysis/plot_all.py --only selection',
    'spineplot/.pixi/envs/default/bin/python3 analysis/plot_all.py --manifest analysis/plot_jobs_v31.toml --only sideband 2>&1 < /dev/null | tee -a plots_v31/logs/x.log',
    'python3 analysis/plotEfficiencyCLI.py --toml a.toml --input f.root --output-dir d --dry-run',
])
def test_automatic(cmd, small_node):
    assert decide(cmd) is None


@pytest.mark.parametrize('cmd', [
    'python3 analysis/plotEfficiencyCLI.py --toml a.toml --input f.root --output-dir d',
    'python3 analysis/plotEfficiency.py',
    'python3 analysis/makeSelectionPlotsCLI.py --toml a.toml --input f.root --output-dir d',
    'python3 analysis/makeNMinus1PlotsCLI.py --toml a.toml --input f.root --output-dir d',
    'python3 analysis/plot_gundam_graphs.py --input g.root',
])
def test_plot_clis_depend_on_memory(cmd, monkeypatch):
    monkeypatch.setattr(guard_bash, 'available_gb', lambda: 344.0)
    assert decide(cmd) is None
    monkeypatch.setattr(guard_bash, 'available_gb', lambda: 9.0)
    d, reason = decide(cmd)
    assert d == 'ask' and 'available' in reason
    monkeypatch.setattr(guard_bash, 'available_gb', lambda: None)   # unknown: ask
    assert decide(cmd)[0] == 'ask'


def test_merge_overwrite_asks(tmp_path):
    out = tmp_path / 'merged_final.root'
    out.write_text('x')
    for cmd in [f'hadd -f {out} a.root b.root', f'~/scripts/mergeFiles.sh -d -l l.txt {out}',
                f"scripts/in_sl7.sh 'hadd {out} a.root'"]:
        d, reason = decide(cmd)
        assert d == 'ask' and 'overwrite' in reason
    # relative paths resolve against the hook's cwd
    assert decide('hadd merged_final.root a.root', cwd=str(tmp_path))[0] == 'ask'
    assert decide('hadd new.root a.root', cwd=str(tmp_path)) is None


# --- read-only commands that only mention guarded names ----------------------

@pytest.mark.parametrize('cmd', [
    'ls -la /pnfs/icarus/scratch/users/micarrig',
    'python3 -m pytest batch -q',
    'git status --short',
    'git ls-remote origin feature/carriganm95_working',
    'jobsub_q --user micarrig',
    'pgrep -fa "[p]lotEfficiencyCLI"',
    'free -g',
    'grep -n medulla CLAUDE.md',
    'cat systematics/src/to_gundam.cc | head',
    'scripts/in_sl7.sh "which root; root-config --version"',
    'cat batch/medulla.py | head -5',
    'grep -rn "to_gundam\\|gundam" plots_v31/logs | head -3',
    'find /exp/icarus/data -maxdepth 3 -iname "*gundam*.root" -printf "%s %p\\n" 2>/dev/null | sort',
    'grep -n hadd CLAUDE.md 2>&1',
    'grep -n jobsub_submit batch/submit.sh',
    'grep pkill scripts/guard_bash.py',
    'git log --oneline -- systematics/add_numi_flux_sys.py',
    'cat scripts/guard_allow.txt',
    "scripts/in_sl7.sh 'grep -n to_gundam CLAUDE.md'",
    "sed -n 590,640p ~/scripts/mergeFiles.sh | grep -n hadd",
    "sed -n '1,20p;40,50p' scripts/guard_bash.py",
    'echo skill',
])
def test_read_only_mentions_pass(cmd):
    assert decide(cmd) is None


@pytest.mark.parametrize('cmd,expected', [
    ('grep x f && rm -rf /tmp/x', 'ask'),                          # each part checked
    (f'ls /tmp; python3 batch/medulla.py -p {P} --launch-jobs', 'ask'),
    (f'ls\npython3 batch/medulla.py -p {P} --launch-jobs', 'ask'),  # newline separates
    (f'python3 analysis/plot_all.py --dry-run && python3 batch/medulla.py -p {P} --launch-jobs', 'ask'),
    ('echo $(kill 1)', 'ask'),                                      # substitution not exempt
    ('echo `kill 1`', 'ask'),
    ('grep x f > out.txt; kill 1', 'ask'),
    ('find /tmp -name "*.root" -exec rm {} +', 'ask'),
    ("sed -n '1p' f; sed -e 's/a/b/' f > g; kill 2", 'ask'),
    (f"scripts/in_sl7.sh 'ls && python3 batch/medulla.py -p {P} -y'", 'ask'),
    ('grep "unbalanced kill 1', 'ask'),                             # unparseable: no exemptions
    (f'echo "unbalanced; python3 batch/medulla.py -p {P}', 'ask'),
])
def test_chains_are_checked_per_part(cmd, expected):
    assert decide(cmd)[0] == expected


# --- owner's allowlist -------------------------------------------------------

ALLOW = [re.compile(r'python3 analysis/makeSelectionPlotsCLI\.py --toml a\.toml'),
         re.compile(r'kill 4242')]


def test_allowlist_skips_ask_rules(small_node):
    assert decide('python3 analysis/makeSelectionPlotsCLI.py --toml a.toml', allow=ALLOW) is None
    assert decide("scripts/in_sl7.sh 'kill 4242'", allow=ALLOW) is None
    assert decide('python3 analysis/makeSelectionPlotsCLI.py --toml a.toml')[0] == 'ask'


@pytest.mark.parametrize('cmd,expected', [
    ('python3 analysis/makeSelectionPlotsCLI.py --toml b.toml', 'ask'),          # not a full match
    ('kill 4242 && rm -rf /tmp/x', 'ask'),                                        # other part
    ('python3 analysis/makeSelectionPlotsCLI.py --toml a.toml; kill 1', 'ask'),
    ('kill 4242 $(rm x)', 'ask'),                                                 # substitution
    (f'kill 4242; rm -rf {P}', 'deny'),                                           # deny still wins
])
def test_allowlist_is_narrow(cmd, expected, small_node):
    assert decide(cmd, allow=ALLOW)[0] == expected


def test_allowlist_file_parsing(tmp_path):
    f = tmp_path / 'allow.txt'
    f.write_text('# comment\n\npython3 x\\.py\n[unclosed\n  # indented comment\n')
    assert [p.pattern for p in guard_bash.load_allowlist(f)] == ['python3 x\\.py']
    assert guard_bash.load_allowlist(tmp_path / 'missing.txt') == []


def test_owner_allowlist_never_covers_deletion():
    for p in guard_bash.load_allowlist():
        assert not p.fullmatch('rm -rf /tmp/x')


# --- the hook as a process -----------------------------------------------------

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
    for cmd in ['ls', 'pkill x', f'python3 batch/medulla.py -p {P} --launch-jobs',
                'python3 batch/medulla.py --dry-run']:
        out = run_hook(json.dumps({'tool_name': 'Bash', 'tool_input': {'command': cmd}}))
        if out:
            assert json.loads(out)['hookSpecificOutput']['permissionDecision'] in ('ask', 'deny')


@pytest.mark.parametrize('payload', ['', 'not json', '{}', '{"tool_name": "Read"}',
                                     '{"tool_name": "Bash", "tool_input": {"command": 5}}'])
def test_hook_fails_open(payload):
    assert run_hook(payload) == ''


def _decision(cmd, mode, cwd=None):
    payload = {'tool_name': 'Bash', 'tool_input': {'command': cmd}}
    if mode is not None:
        payload['permission_mode'] = mode
    if cwd is not None:
        payload['cwd'] = cwd
    out = run_hook(json.dumps(payload))
    return json.loads(out)['hookSpecificOutput']['permissionDecision'] if out else None


@pytest.mark.parametrize('mode', ['auto', 'bypassPermissions', 'dontAsk'])
def test_ask_becomes_deny_when_nobody_is_prompted(mode):
    for cmd in ['kill 123', f'python3 batch/medulla.py -p {P} -e icarus --launch-jobs --tag abc',
                f'python3 batch/medulla.py -p {P} -e icarus -Z -y',
                'python3 analysis/plot_all.py --only nminus1; kill 123']:
        assert _decision(cmd, mode) == 'deny'
    # Dry-runs, status checks, local stages and ordinary commands are untouched.
    for cmd in [f'python3 batch/medulla.py -p {P} -e icarus -y --launch-jobs --dry-run',
                f'python3 batch/medulla.py -p {P} -e icarus -Z < /dev/null',
                'python3 analysis/plot_all.py --only nminus1',
                'python3 batch/merge_list.py P --output m.txt',
                'build/systematics/to_gundam c.toml m.root', 'ls']:
        assert _decision(cmd, mode) is None


def test_hook_uses_payload_cwd(tmp_path):
    (tmp_path / 'exists.root').write_text('x')
    assert _decision('hadd exists.root a.root', None, cwd=str(tmp_path)) == 'ask'
    assert _decision('hadd new.root a.root', None, cwd=str(tmp_path)) is None


@pytest.mark.parametrize('mode', [None, 'default', 'acceptEdits', 'plan'])
def test_ask_stays_ask_when_owner_is_prompted(mode):
    assert _decision(f'python3 batch/medulla.py -p {P} -e icarus --launch-jobs', mode) == 'ask'
    assert _decision('kill 123', mode) == 'ask'
