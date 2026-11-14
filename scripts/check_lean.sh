#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    printf 'Usage: %s [--chain] [FILE ...]\n' "$0"
    printf 'Build the pinned Lean project, then check any supplied Lean source files.\n'
    printf 'Both .lean and .txt files are accepted. With no files, check the project only.\n'
    printf 'With --chain the files are one chain of modules, checked in order by one Lean process\n'
    printf 'each: a file may import the files before it by name (NAME.lean as NAME). Each but the\n'
    printf 'last is compiled, under the same flags, to an .olean in a temporary directory that only\n'
    printf 'the later files of this run read.\n'
    exit 0
fi
chain=0
if [[ "${1:-}" == "--chain" ]]; then
    chain=1
    shift
fi
if [[ "${1:-}" == "--" ]]; then
    shift
fi

repo_dir="$(CDPATH= cd -- "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
project_dir="$repo_dir/lean"
toolchain="$(< "$project_dir/lean-toolchain")"

if command -v elan >/dev/null 2>&1; then
    elan_command="$(command -v elan)"
elif [[ -x "${HOME:-}/.elan/bin/elan" ]]; then
    elan_command="${HOME}/.elan/bin/elan"
else
    printf 'check_lean.sh: elan was not found. Install elan; see lean/README.md.\n' >&2
    exit 127
fi

for file in "$@"; do
    if [[ ! -f "$file" || ! -r "$file" ]]; then
        printf 'check_lean.sh: cannot read file: %s\n' "$file" >&2
        exit 2
    fi
done
if (( chain )); then
    # An imported file's module name is its file name: an identifier, so that `import NAME` finds it.
    for file in "${@:1:$(( $# > 0 ? $# - 1 : 0 ))}"; do
        if [[ ! "$(basename -- "$file")" =~ ^[A-Za-z_][A-Za-z0-9_]*\.lean$ ]]; then
            printf 'check_lean.sh: a chained file that later files import must be named IDENTIFIER.lean: %s\n' "$file" >&2
            exit 2
        fi
    done
fi

caller_dir="$PWD"
cd -- "$project_dir"
# Z3_LEAN_SKIP_BUILD=1: the caller has already built the project in this checkout (the Python
# tests set it after their first check).
if [[ "${Z3_LEAN_SKIP_BUILD:-}" != "1" ]]; then
    "$elan_command" run "$toolchain" lake build
fi
if (( chain )); then
    # The .olean files of the chain: written by this run's Lean processes, read by its later ones,
    # and removed with the run. `lake env` keeps LEAN_PATH after the project's own search path.
    modules="$(mktemp -d)"
    trap 'rm -rf -- "$modules"' EXIT
    export LEAN_PATH="$modules${LEAN_PATH:+:$LEAN_PATH}"
fi
remaining=$#
for file in "$@"; do
    case "$file" in
        /*) ;;
        *) file="$caller_dir/$file" ;;
    esac
    remaining=$(( remaining - 1 ))
    printf 'Checking %s\n' "$file"
    if (( chain && remaining )); then
        name="$(basename -- "$file")"
        "$elan_command" run "$toolchain" lake env lean --trust=0 -DwarningAsError=true --memory=8192 \
            --root="$(dirname -- "$file")" -o "$modules/${name%.lean}.olean" "$file"
    else
        "$elan_command" run "$toolchain" lake env lean --trust=0 -DwarningAsError=true --memory=8192 "$file"
    fi
done
printf 'Lean checks passed.\n'
