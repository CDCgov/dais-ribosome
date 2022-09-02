# shellcheck shell=bash
shopt -s expand_aliases

PROG=${IFX_PROG:-$(basename "$0")}
DEBUG=${IFX_DEBUG:-0}

function die() {
    local t=$(date +"%Y-%m-%d %k:%M:%S")
    local pad=" $(printf "%*s" "$(((SHLVL - 2) * 2))" "")"

    echo -e "[$t]${pad}$PROG ERROR :: $1" 1>&2
    exit 1
}

function warn() {
    local t=$(date +"%Y-%m-%d %k:%M:%S")
    local pad=" $(printf "%*s" "$(((SHLVL - 2) * 2))" "")"
    echo -e "[$t]${pad}$PROG WARNING :: $1" 1>&2
    return 1
}

alias bail='return 1'

function time_stamp() {
    [ "${DEBUG:-0}" -eq "1" ] \
        && echo -e "  FUNCTION TRACE:" "${FUNCNAME[@]}"

    local t=$(date +"%Y-%m-%d %k:%M:%S")
    local pad=" $(printf "%*s" "$(((SHLVL - 2) * 2))" "")"

    echo -e "[$t]${pad}$PROG :: $1"
}

function finish() {
    time_stamp "$1"
    [ "$(type -t unlock)" == "function" ]   && unlock
    [ "$(type -t cleanup)" == "function" ]  && cleanup
    exit 0
}

function time_stamp_V() {
    [ "${DEBUG:-0}" -eq "1" ] \
        && echo -e "  FUNCTION TRACE:" "${FUNCNAME[@]}"

    if [ "${VERBOSE:-0}" -eq "1" ]; then
        time_stamp "  $1"
    fi
}

# Idea courtesy: stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
function check_program() {
    command -v "$1" > /dev/null 2>&1 || warn "Program '$1' not found, add to the PATH or install."
}
