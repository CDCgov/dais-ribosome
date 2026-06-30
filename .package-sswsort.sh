#!/bin/bash

default_sswsort=v2.2.4

version=${PINNED_SSWSORT:-$default_sswsort}
folder=sswsort-${version#v}
archive=${folder}.zip
url=https://github.com/CDCgov/sswsort/archive/refs/tags/${version}.zip

curl -L --proto '=https' --tlsv1.2 -sSf --output "$archive" "$url" \
    && unzip -o "$archive" \
    && rm "$archive" \
    && mv "${folder}/sswsort_res" . \
    && rm -r "${folder}"
