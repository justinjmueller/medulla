#!/bin/bash

## Container start up (for reference; don't run here)
# sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh

function set_samweb() {
    export SAM_EXPERIMENT=$(id -ng) 
    export SAM_GROUP=$(id -ng)
    export SAM_STATION=$(id -ng) 
    export SAM_WEB_HOST=sam$(id -ng).fnal.gov 
    export IFDH_BASE_URI=http://sam$(id -ng).fnal.gov:8480/sam/$exp/api/
    setup sam_web_client
}

function get_bearer_token() {
    export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u)
    htgettoken -a htvaultprod.fnal.gov --vaulttokenttl=1d --vaulttokenminttl=12h -i $(id -ng)
}

function apply_tweaks() {
    # Fix the terminal type in screen sessions. This is important if
    # you are using a terminal multiplexer like screen or tmux!
    export TERM=xterm

    # Enable color support for ls and grep if available
    if command -v dircolors >/dev/null 2>&1; then
        # Set up LS_COLORS
        eval "$(dircolors -b)"
        alias ls='ls --color=auto'
        alias grep='grep --color=auto'
        alias egrep='egrep --color=auto'
        alias fgrep='fgrep --color=auto'
    fi

    # Optional: colorized prompt (bash only)
    if [[ $- == *i* ]] && [[ -n "$PS1" ]]; then
        # Set a colored prompt: user@host:cwd$
        #PS1='\[\e[0;32m\]\u@\h\[\e[0m\]:\[\e[0;34m\]\w\[\e[0m\]\$ '
        PS1='\[\e[0;32m\]$>\[\e[0m\] '
    fi
}

# Setup CVMFS area
source /cvmfs/$(id -ng).opensciencegrid.org/products/$(id -ng)/setup_$(id -ng).sh

# Setup the required dependencies
setup sbnana v10_01_04 -q e26:prof
setup cmake v3_27_4

# Setup SAMWeb and the proxy certificate
set_samweb
get_bearer_token

# Apply tweaks to the environment
apply_tweaks
