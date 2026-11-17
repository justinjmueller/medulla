"""
Data, rock and detector-variation samples get no systematics output: merge B
and phase 1 read their output_jobid*.root, phase 2 only takes 'nominal', and
the detector-variation CAFs carry no rec.mc.nu.wgt (add_weights aborted).

- create_systematics_cfg writes no [[tree]] blocks for them;
- submit.sh exits before run_systematics for them;
- check_project_status does not flag their leftover systematics files (stubs
  or recovered files from older projects), so a -y run cannot delete the good
  selection output with them.
"""
import subprocess
from pathlib import Path

import pytest

from test_dry_run import GOOD, STUB, make_project, job_statuses, work, reset_flags  # noqa: F401
import utilities

BATCH_DIR = Path(__file__).parent

TREE = {
    'name': 'selected',
    'add_systematics': True,
    'branch': [{'name': 'neutrino_id', 'type': 'true'},
               {'name': 'neutrino_energy', 'type': 'mctruth'}],
}


def base_cfg():
    return {'input': {'use_additional_hash': True}, 'output': {}}


def test_cfg_skips_no_systematics_tags():
    samples = [
        {'name': 'NuMIFull', 'ismc': True, 'tag': 'nominal'},
        {'name': 'data', 'ismc': False, 'tag': 'data'},
        {'name': 'coherent_high', 'ismc': True, 'tag': 'detector_variation'},
        {'name': 'cv', 'ismc': True, 'tag': 'detector_variation'},
        {'name': 'rock', 'ismc': True, 'tag': 'rock'},
        {'name': 'NuMI_nue', 'ismc': True, 'tag': 'nue_enhanced'},
        {'name': 'g4', 'ismc': True, 'tag': 'g4'},
    ]
    cfg = utilities.create_systematics_cfg(base_cfg(), [TREE], samples)
    origins = {t['origin']: t['action'] for t in cfg['tree']}
    assert origins == {
        'events/NuMIFull/selected': 'add_weights',
        'events/NuMI_nue/selected': 'add_weights',
        'events/g4/selected': 'add_weights',
    }


def test_status_ignores_no_systematics_sys_files(tmp_path, work):
    project = make_project(
        tmp_path,
        jobs=[(0, 'detector_variation', 'submitted'), (1, 'nominal', 'submitted'),
              (2, 'data', 'submitted'), (3, 'rock', 'submitted')],
        outputs={
            'output_jobid0000.root': GOOD,
            'output_systematics_jobid0000.root': STUB,   # aborted run_systematics: ignore
            'output_jobid0001.root': GOOD,
            'output_systematics_jobid0001.root': STUB,   # nominal: still a real failure
            'output_jobid0002.root': GOOD,
            'output_systematics_jobid0002.root': STUB,   # data: ignore
            'output_jobid0003.root': GOOD,               # rock: no systematics file at all
        },
    )
    utilities.set_assume_yes(True)
    utilities.check_project_status(project)
    out = project / 'output'
    assert (out / 'output_jobid0000.root').exists()
    assert (out / 'output_systematics_jobid0000.root').exists()
    assert job_statuses(project)[0] == 'completed'
    assert not (out / 'output_jobid0001.root').exists()
    assert not (out / 'output_systematics_jobid0001.root').exists()
    assert job_statuses(project)[1] == 'pending'
    assert (out / 'output_jobid0002.root').exists()
    assert job_statuses(project)[2] == 'completed'
    assert job_statuses(project)[3] == 'completed'


@pytest.mark.parametrize('tag,runs', [
    ('detector_variation', False), ('data', False), ('rock', False),
    ('nominal', True), ('nue_enhanced', True), ('g4', True),
])
def test_submit_sh_skips_run_systematics(tmp_path, tag, runs):
    """Run only submit.sh's systematics block, with run_systematics and ifdh stubbed."""
    script = (BATCH_DIR / 'submit.sh').read_text()
    start = script.index('# Data, rock and detector-variation samples get no')
    block = script[start:]
    (tmp_path / 'job_config.toml').write_text(f'[[sample]]\nname = "s"\ntag = "{tag}"\nismc = true\n')
    (tmp_path / 'systematics').mkdir()
    marker = tmp_path / 'ran'
    rs = tmp_path / 'systematics' / 'run_systematics'
    rs.write_text(f'#!/bin/bash\ntouch {marker}\n')
    rs.chmod(0o755)
    prelude = 'ifdh() { :; }\nJOBID=7\nPROJECT=/nonexistent\n'
    res = subprocess.run(['bash', '-c', prelude + block], cwd=tmp_path, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert marker.exists() == runs
