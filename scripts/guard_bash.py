#!/usr/bin/env python3
"""
Claude Code PreToolUse hook for the Bash tool: enforce the rules in CLAUDE.md.

Reads the hook JSON on stdin ({"tool_name": "Bash", "tool_input": {"command": ...}})
and, for commands that match a rule, prints a PreToolUse decision:

  deny : never allowed from an agent (kills the shell, hidden htgettoken, raw
         jobsub_submit, deleting project outputs)
  ask  : allowed only if the owner approves (anything that submits, deletes,
         or loads large inputs, unless it is a --dry-run / --help)

In auto mode (and the other modes that approve without prompting) an "ask"
would be approved without the owner seeing it, so there it becomes "deny": the
owner switches out of auto mode (shift+tab) to approve the command. The one
exception is a plotting run in auto mode: it stays "ask" and the auto-mode
classifier allows it only if the owner's latest message approves that run
(autoMode.soft_deny in .claude/settings.local.json).

It never prints "allow", so it can only make the normal permission flow
stricter. If it cannot parse its input it prints nothing (fail open), so a bug
here never blocks ordinary commands.

Tested by batch/test_guard_bash.py.
"""
import json
import re
import sys

# --- deny ------------------------------------------------------------------

DENY = [
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

# --- ask -------------------------------------------------------------------

HARMLESS_FLAGS = re.compile(r'(--dry-run|--help|(?<![\w-])-h(?![\w-]))')

# (pattern, reason). Each is skipped when the command also has --dry-run/--help/-h.
ASK_UNLESS_DRY = [
    (r'batch/medulla\.py|(?<![\w/])medulla\.py',
     'medulla.py can create projects, delete outputs and submit jobs. Run it with --dry-run first and show the owner.'),
    (r'plotEfficiency(CLI)?\.py|make\w*Plots(CLI)?\.py|plot_gundam_graphs\.py|plot_all\.py',
     'Plotting scripts load 30 GB+ inputs into memory. Use --dry-run (if available) and ask first (CLAUDE.md rule 10).'),
]
# The plotting rule above; in auto mode it is left to the auto-mode classifier
# (autoMode.soft_deny in .claude/settings.local.json: the owner's latest message
# must approve the run) instead of being blocked, so the owner can approve from
# a phone. Everything else asked about is still blocked in auto mode.
PLOTTING = ASK_UNLESS_DRY[1][0]
ASK_UNLESS_DRY += [
    (r'batch/merge_list\.py|(?<![\w/])merge_list\.py',
     'merge_list.py copies the whole project.db (hundreds of MB) and stats every output on PNFS. Run it with --dry-run first and show the owner.'),
]

# Heavy or output-writing; no dry-run mode exists, so always ask.
ASK_ALWAYS = [
    (r'mergeFiles\.sh|(?<![\w/])hadd\b',
     'Merging reads and writes very large ROOT files. Ask first.'),
    (r'add_numi_flux_sys\.py',
     'add_numi_flux_sys.py edits the merged file in place and is slow. Ask first.'),
    (r'(?<!\w)to_gundam(?![\w.])',
     'to_gundam writes output.root in the current directory. Ask first.'),
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


# Modes in which an "ask" is approved without a prompt reaching the owner.
NO_PROMPT_MODES = {'auto', 'bypassPermissions', 'dontAsk'}


def decide(command: str):
    """Return (decision, reason) or None."""
    for pat, reason in DENY:
        if re.search(pat, command):
            return 'deny', reason
    harmless = bool(HARMLESS_FLAGS.search(command))
    for pat, reason in ASK_UNLESS_DRY:
        if re.search(pat, command) and not harmless:
            return 'ask', reason
    for pat, reason in ASK_ALWAYS:
        if re.search(pat, command):
            return 'ask', reason
    return None


def only_plotting(command: str) -> bool:
    """True if the plotting rule is the only ask rule the command matches."""
    if not re.search(PLOTTING, command):
        return False
    others = [p for p, _ in ASK_UNLESS_DRY if p != PLOTTING] + [p for p, _ in ASK_ALWAYS]
    return not any(re.search(p, command) for p in others)


def main():
    try:
        data = json.load(sys.stdin)
        if data.get('tool_name') != 'Bash':
            return
        command = data.get('tool_input', {}).get('command', '')
        if not isinstance(command, str):
            return
    except Exception:
        return
    result = decide(command)
    if result is None:
        return
    decision, reason = result
    mode = data.get('permission_mode')
    plotting_in_auto = mode == 'auto' and only_plotting(command)
    if decision == 'ask' and mode in NO_PROMPT_MODES and not plotting_in_auto:
        decision = 'deny'
        reason += (f" Blocked because the session is in {data['permission_mode']} mode, which would "
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
