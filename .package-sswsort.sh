#!/bin/bash

default_sswsort=v2.0.0

version=${PINNED_SSWSORT:-$default_sswsort}
folder=sswsort-${version#v}
archive=${folder}.zip
url=https://github.com/CDCgov/sswsort/archive/refs/tags/${version}.zip

wget "$url" --output-document "$archive" \
    && unzip -o "$archive" \
    && rm "$archive" \
    && mv "${folder}/sswsort_res" . \
    && rm -r "${folder}"
