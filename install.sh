#!/bin/bash

app=ribosome

if command -v "cargo" > /dev/null 2>&1; then
    if [ -n "$1" ]; then
        if grep -Eq "not a recognized processor" <(rustc -C "target-cpu=$1" --print host-tuple 2>&1); then
            echo "WARNING: Cannot find target-cpu '$1', skipping"
        else
            export RUSTFLAGS="${RUSTFLAGS:+$RUSTFLAGS }-C target-cpu=$1"
        fi
    fi

    cargo +nightly build --profile prod \
        && cp target/prod/$app . \
        && ./.package-sswsort.sh \
        || {
            echo "ERROR: Installation failed!"
            exit 1
        }
else
    echo "WARNING: Rust / cargo not found!"
    read -r -p "INSTALLER: Download the ribosome binary instead? [Y/n]: " yn
    yn=${yn:-y}
    case $yn in
        [Yy]*)
            ./.download-ribosome-binary.sh
            ;;
        [Nn]*)
            echo "INSTALLER: Exiting."
            exit 0
            ;;
        *)
            echo "ERROR: Invalid input. Exiting."
            exit 1
            ;;
    esac
fi
