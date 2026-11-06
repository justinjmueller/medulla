#!/usr/bin/env python3
"""
Claude Code PreToolUse hook for the Bash tool: enforce the rules in CLAUDE.md.

Reads the hook JSON on stdin ({"tool_name": "Bash", "tool_input": {"command": ...},
"cwd": ..., "permission_mode": ...}) and, for commands that match a rule, prints
a PreToolUse decision:

  deny : never allowed from an agent (kills the shell, hidden htgettoken, raw
         jobsub_submit, deleting project outputs, editing guard_allow.txt)
  ask  : the owner approves it at a real prompt. Only for: submitting grid
         jobs with medulla.py, medulla.py --yes (deletes outputs unprompted),
         plotting CLIs on a node with too little memory, overwriting an
         existing merge output, add_numi_flux_sys.py without --output,
         deleting/moving files, kill, htgettoken.

Everything else runs without asking (owner's decision, 2026-10-23): every
--dry-run/--help, the medulla.py status check, project creation (-c),
interactive phase 1, a single --test-job (with -y), merge_list.py, merges into a new file, stage 5 with
--output, to_gundam, plot_all.py (it has its own memory gate) and the plotting
CLIs on a node with at least PLOT_MIN_AVAIL_GB available.

In auto mode (and the other modes that approve without prompting) an "ask"
would be approved without the owner seeing it, so there it becomes "deny": the
owner switches out of auto mode (shift+tab) to approve the command.

The command is split into simple commands (at unquoted ; & | && || and
newlines) and each part is checked on its own, so a --dry-run in one part does
not excuse another part. A part is exempt from every rule when it is a
read-only command (grep, ls, cat, sed -n '<range>p', find without -exec/-delete,
...) with no output redirection or command substitution: mentioning a guarded
name is not running it. `scripts/in_sl7.sh '<cmd>'` is checked as `<cmd>`.

scripts/guard_allow.txt is the owner's list of commands marked as not causing
memory issues: one Python regex per line, full-matched against a part. A part
that matches skips the ask rules (never the deny rules). Only the owner edits
that file; from Bash any write to it is denied, and settings.local.json denies
Edit/Write on it. If the quoting cannot be parsed, the whole command is checked
as one part with no exemptions.

It never prints "allow", so it can only make the normal permission flow
stricter. If it cannot parse its input it prints nothing (fail open), so a bug
here never blocks ordinary commands.

Tested by batch/test_guard_bash.py.
"""
import json
import os
import re
import shlex
import sys
from pathlib import Path

ALLOW_FILE = Path(__file__).resolve().parent / 'guard_allow.txt'

# Plotting CLIs other than plot_all.py load whole trees (30 GB+ inputs) and have
# no memory gate of their own; below this much available memory they ask.
PLOT_MIN_AVAIL_GB = 64

# --- deny ------------------------------------------------------------------

DENY = [
    (r'guard_allow\.txt',
     'Only the owner edits scripts/guard_allow.txt. Propose the line to the owner instead.'),
    (r'(?<![\w./-])(pkill|killall)\b',
     'pkill/killall by name can match the agent\'s own shell. Use pgrep -fa "[n]ame" and kill by PID.'),
    (r'(?<![\w./-])(source|\.)\s+\S*setup_icarus_env\.sh',
     'setup_icarus_env.sh runs htgettoken and sets "set -euo pipefail". Use scripts/in_sl7.sh for the environment.'),
    (r'setup_icarus_env\.sh',
     'setup_icarus_env.sh runs htgettoken. Use scripts/in_sl7.sh for the environment.'),
    (r'(?<![\w./-])jobsub_(submit|rm|hold|release)\b',
     'Raw jobsub commands are not allowed. Submit through batch/medulla.py (show its --dry-run output first).'),
    (r'(?<![\w./-])(rm|unlink|rmdir|mv|shred|truncate)\b[^;&|\n]*(/pnfs/|nueCCInclusiveSpine)',
     'Deleting or moving files in /pnfs or a project directory needs the owner to do it (CLAUDE.md rule 2/9).'),
    (r'(?<![\w./-])find\b[^;&|\n]*(/pnfs/|nueCCInclusiveSpine)[^;&|\n]*(-delete|-exec(dir)?\s+(rm|mv|unlink|shred)\b)',
     'find -delete / -exec rm on /pnfs or a project directory needs the owner to do it (CLAUDE.md rule 2/9).'),
]

# --- ask (pattern rules) ---------------------------------------------------

HARMLESS_FLAGS = re.compile(r'(--dry-run|--help|(?<![\w-])-h(?![\w-]))')

ASK_ALWAYS = [
    (r'(?<![\w./-])(rm|mv|unlink|rmdir|shred|truncate)\b',
     'Deleting, moving or truncating files needs approval (CLAUDE.md rule 9). Name the exact files.'),
    (r'(?<![\w./-])find\b[^;&|\n]*(-delete|-exec(dir)?\s+(rm|mv|unlink|shred)\b)',
     'find -delete / -exec rm deletes files. Needs approval (CLAUDE.md rule 9).'),
    (r'rmtree|os\.(remove|unlink|rmdir)\b|\.unlink\(|shutil\.move|os\.(rename|replace)\b',
     'Deleting or moving files from inline Python needs approval (CLAUDE.md rule 9).'),
    (r'(?<![\w./-])kill\b',
     'Killing a process (even by PID) can stop a run the owner is waiting on. Ask first.'),
    (r'(?<![\w./-])htgettoken\b',
     'Token refresh may need a browser and is the owner\'s call (CLAUDE.md: ask before htgettoken).'),
]

# --- ask (rules that look at the arguments) ---------------------------------

MEDULLA = re.compile(r'(^|/)medulla\.py$')
PLOT_CLI = re.compile(r'(^|/)(plotEfficiency(CLI)?|make\w*Plots(CLI)?|plot_gundam_graphs)\.py$')
SUBMIT_FLAGS = {'--launch-jobs', '--test-job', '--variation-phase2', '--relaunch-missing'}

REASON_SUBMIT = ('medulla.py would submit grid jobs. Show the owner the --dry-run output '
                 '(project, jobs, tag/SHA, memory/disk, jobsub_submit command) and get a yes.')
REASON_YES = ('medulla.py --yes deletes stub/zombie outputs without prompting (CLAUDE.md rule 2). '
              'Run the read-only status check first and get the owner\'s yes for that run.')
REASON_PLOT_MEM = ('This plotting CLI loads whole trees and this node has only {:.0f} GB available '
                   '(needs {} GB). Run it on the big-memory node, or use plot_all.py (memory-gated).')
REASON_OVERWRITE = 'The merge would overwrite the existing file {}. Ask the owner first.'
REASON_FLUX = ('add_numi_flux_sys.py without --output edits the merged file in place (that is how v31\'s '
               'first run corrupted it). Pass --output <new>.root.')


def available_gb():
    """MemAvailable from /proc/meminfo in GB; None if it cannot be read."""
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) / 1024 ** 2
    except (OSError, ValueError, IndexError):
        pass
    return None


def _script(words):
    """The script or binary a command runs: skips python/bash/sh/env wrappers."""
    for w in words:
        base = os.path.basename(w)
        if w.startswith('-') or re.match(r'^(python[\d.]*|bash|sh|env)$', base):
            continue
        return w
    return ''


def _short_flag(words, letter):
    """True if a single-dash flag cluster (-y, -Zy) contains `letter`."""
    return any(re.fullmatch(r'-[A-Za-z]+', w) and letter in w[1:] for w in words)


def _exists(path, cwd):
    p = Path(os.path.expanduser(path))
    if not p.is_absolute():
        p = Path(cwd or os.getcwd()) / p
    return p.exists()


def _arg_hits(seg, words, cwd):
    """[(pattern_id, reason)] from the rules that need the arguments."""
    if not words:
        return []
    harmless = bool(HARMLESS_FLAGS.search(seg))
    script = _script(words)
    name = os.path.basename(script)
    hits = []
    if MEDULLA.search(script) and not harmless:
        interactive_p1 = '--variation-phase1' in words and '--variation-interactive' in words
        grid_p1 = '--variation-phase1' in words and not interactive_p1
        # --test-job/-T submits exactly one job (selection, or phase 2 when combined
        # with --variation-phase2); the owner allows it without approval (2026-10-23).
        test_job = (('--test-job' in words or _short_flag(words, 'T'))
                    and '--launch-jobs' not in words and not grid_p1)
        submit_flags = SUBMIT_FLAGS - {'--test-job', '--variation-phase2', '--relaunch-missing'} \
            if test_job else SUBMIT_FLAGS
        if submit_flags & set(words) or grid_p1 or (_short_flag(words, 'r') and not test_job):
            hits.append(('medulla-submit', REASON_SUBMIT))
        if ('--yes' in words or _short_flag(words, 'y')) and not (interactive_p1 or test_job):
            hits.append(('medulla-yes', REASON_YES))
    if PLOT_CLI.search(script) and not harmless:
        avail = available_gb()
        if avail is None or avail < PLOT_MIN_AVAIL_GB:
            hits.append(('plot-mem', REASON_PLOT_MEM.format(avail or 0, PLOT_MIN_AVAIL_GB)))
    if name in ('hadd', 'mergeFiles.sh'):
        roots = [w for w in words[1:] if not w.startswith('-') and w.endswith('.root')]
        if roots:
            out = roots[0] if name == 'hadd' else roots[-1]
            if _exists(out, cwd):
                hits.append(('merge-overwrite', REASON_OVERWRITE.format(out)))
    if name == 'add_numi_flux_sys.py' and not harmless and not any(
            w == '--output' or w.startswith('--output=') for w in words):
        hits.append(('flux-inplace', REASON_FLUX))
    return hits


# Modes in which an "ask" is approved without a prompt reaching the owner.
NO_PROMPT_MODES = {'auto', 'bypassPermissions', 'dontAsk'}


# --- splitting and exemptions -----------------------------------------------

# Commands that only look. Exempt only without output redirection or substitution.
READ_ONLY = {
    'grep', 'egrep', 'fgrep', 'rg', 'ls', 'cat', 'head', 'tail', 'wc', 'stat', 'file',
    'du', 'df', 'echo', 'printf', 'pgrep', 'which', 'type', 'free', 'uptime', 'pwd',
    'basename', 'dirname', 'realpath', 'readlink', 'diff', 'cmp', 'md5sum', 'sha256sum',
    'cut', 'tr', 'column', 'tree', 'id', 'hostname', 'date', 'true', 'cd', 'find',
}
READ_ONLY_GIT = {'status', 'log', 'diff', 'show', 'grep', 'ls-remote', 'rev-parse', 'ls-files', 'blame'}
FIND_ACTIONS = re.compile(r'(?<![\w-])-(delete|exec|execdir|ok|okdir|fprint|fprint0|fprintf|fls)\b')
SED_PRINT = re.compile(r'(\d+|\$)(,(\d+|\$))?p(;(\d+|\$)(,(\d+|\$))?p)*')
SUBSTITUTION = re.compile(r'\$\(|`|[<>]\(')
HARMLESS_REDIRECT = re.compile(r'\d*>>?\s*(/dev/null\b|&\d)|\d*<\s*/dev/null\b')


def split_segments(command: str):
    """Split at unquoted ; & | and newlines. None if the quoting is unbalanced."""
    segs, cur, quote, i, n = [], [], None, 0, len(command)
    while i < n:
        c = command[i]
        if quote:
            cur.append(c)
            if c == '\\' and quote == '"' and i + 1 < n:
                cur.append(command[i + 1])
                i += 1
            elif c == quote:
                quote = None
        elif c == '\\' and i + 1 < n:
            cur += [c, command[i + 1]]
            i += 1
        elif c in '\'"':
            quote = c
            cur.append(c)
        elif c == '&' and ((cur and cur[-1] in '<>') or command[i + 1:i + 2] == '>'):
            cur.append(c)  # part of a redirection: 2>&1, &>
        elif c in ';&|\n':
            segs.append(''.join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    if quote:
        return None
    segs.append(''.join(cur))
    return [s.strip() for s in segs if s.strip()]


def _words(seg: str):
    try:
        words = shlex.split(seg)
    except ValueError:
        return None
    while words and re.match(r'^\w+=', words[0]):  # FOO=bar cmd
        words.pop(0)
    return words


def is_read_only(seg: str) -> bool:
    if SUBSTITUTION.search(seg) or '>' in HARMLESS_REDIRECT.sub('', seg):
        return False
    words = _words(seg)
    if not words:
        return False
    name = words[0]
    if name == 'git':
        return len(words) > 1 and words[1] in READ_ONLY_GIT
    if name == 'find':
        return not FIND_ACTIONS.search(seg)
    if name == 'sed':  # only the print-a-range form: sed -n '10,20p' file...
        return len(words) >= 3 and words[1] == '-n' and bool(SED_PRINT.fullmatch(words[2]))
    return name in READ_ONLY


def load_allowlist(path=None):
    """The owner's regexes from guard_allow.txt; bad or missing lines are ignored."""
    pats = []
    try:
        lines = Path(path or ALLOW_FILE).read_text().splitlines()
    except OSError:
        return pats
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            pats.append(re.compile(line))
        except re.error:
            pass
    return pats


def is_allowlisted(seg: str, allow) -> bool:
    if SUBSTITUTION.search(seg):
        return False
    return any(p.fullmatch(seg) for p in allow)


def _segment_hits(seg: str, allow, cwd, depth=0):
    """[(decision, pattern, reason)] for one simple command."""
    words = _words(seg)
    if (depth < 3 and words and len(words) == 2 and words[0].endswith('in_sl7.sh')
            and not SUBSTITUTION.search(seg)):
        return _command_hits(words[1], allow, cwd, depth + 1)
    if is_read_only(seg):
        return []
    hits = [('deny', p, r) for p, r in DENY if re.search(p, seg)]
    if hits or is_allowlisted(seg, allow):
        return hits
    hits += [('ask', p, r) for p, r in _arg_hits(seg, words, cwd)]
    hits += [('ask', p, r) for p, r in ASK_ALWAYS if re.search(p, seg)]
    return hits


def _command_hits(command: str, allow, cwd, depth=0):
    segs = split_segments(command)
    if segs is None:  # unparseable: whole command, no exemptions
        hits = [('deny', p, r) for p, r in DENY if re.search(p, command)]
        hits += [('ask', p, r) for p, r in ASK_ALWAYS if re.search(p, command)]
        if not hits and re.search(r'medulla\.py|plot\w*\.py|Plots\w*\.py|hadd|mergeFiles\.sh|'
                                  r'add_numi_flux_sys\.py', command):
            hits.append(('ask', 'unparsed', 'Could not parse the quoting of this command; ask first.'))
        return hits
    return [h for s in segs for h in _segment_hits(s, allow, cwd, depth)]


def decide(command: str, allow=None, cwd=None):
    """Return (decision, reason) or None."""
    hits = _command_hits(command, load_allowlist() if allow is None else allow, cwd)
    for d in ('deny', 'ask'):
        for decision, _, reason in hits:
            if decision == d:
                return decision, reason
    return None


def main():
    try:
        data = json.load(sys.stdin)
        if data.get('tool_name') != 'Bash':
            return
        command = data.get('tool_input', {}).get('command', '')
        if not isinstance(command, str):
            return
        cwd = data.get('cwd') if isinstance(data.get('cwd'), str) else None
    except Exception:
        return
    result = decide(command, cwd=cwd)
    if result is None:
        return
    decision, reason = result
    mode = data.get('permission_mode')
    if decision == 'ask' and mode in NO_PROMPT_MODES:
        decision = 'deny'
        reason += (f" Blocked because the session is in {mode} mode, which would "
                   "approve it without asking. Ask the owner; they switch out of that mode "
                   "(shift+tab) to approve it.")
    print(json.dumps({
        'hookSpecificOutput': {
            'hookEventName': 'PreToolUse',
            'permissionDecision': decision,
            'permissionDecisionReason': f'[guard_bash] {reason}',
        }
    }))


if __name__ == '__main__':
    main()
