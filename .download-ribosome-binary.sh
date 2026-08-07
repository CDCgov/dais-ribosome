#!/bin/bash

app=dais-ribosome
version=$(git describe --tags --abbrev=0) || exit 1
archive=ribosome-cli-download.tar.gz

if [[ "$(uname -o)" == "Darwin" ]]; then
    os_arch="apple-universal"
elif [[ "$(uname -m)" == "x86_64" ]]; then
    os_arch="linux-x86_64"
elif [[ "$(uname -m)" == "aarch64" ]]; then
    os_arch="linux-aarch64"
else
    echo "ERROR: unsupported OS / architecture"
    exit 1
fi

if [ -d "./sswsort_res" ]; then
    echo "INSTALLER: Removing previous './sswsort_res/'"
    rm -rf ./sswsort_res
fi

url=https://github.com/CDCgov/${app}/releases/download/${version}/${app}-cli-${os_arch}-${version}.tar.gz
echo "INSTALLER: Downloading '$url'"
curl -L --proto '=https' --tlsv1.2 -sSf --output "$archive" "$url" \
    && tar xzf $archive \
    && mv "${app}-cli-${version}"/{sswsort_res,ribosome} . \
    && rm -rf "$archive" "${app}-cli-${version}" \
    && echo "INSTALLER: Download success!" \
    || {
        echo "ERROR: Download failed. You may require curl or there was a network issue."
        exit 1
    }
