#!/bin/bash

VERSION="v1.0"
COMP="clean-up"

export IFX_LOCAL_PROCS=48

if [ "$#" == "2" ]; then
    flu_data="cat $1"
    cov_data="cat $2"
else
    flu_data="lib/editMSA/ordinalHeaders.pl spec/INFLUENZA.refs"
    cov_data="lib/editMSA/ordinalHeaders.pl spec/BETACORONAVIRUS.refs"
fi

for tag in "$VERSION" "$COMP"; do
    git checkout $tag > /dev/null \
        || {
            echo "Git failure could be do to unstashed or uncommited work. Stash/commit and then try again."
            exit 1
        }
    echo "Running '$tag'" \
        && time ./ribosome --module INFLUENZA <($flu_data) \
            test-flu-$tag.seq test-flu-$tag.ins test-flu-$tag.del test-flu-$tag.gen \
        && time ./ribosome --module BETACORONAVIRUS <($cov_data) \
            test-cov-$tag.seq test-cov-$tag.ins test-cov-$tag.del test-cov-$tag.gen

done

for mod in flu cov; do
    for ext in seq ins del gen gen.ins gen.del; do
        echo -n "Checking $mod / $ext" \
            && cmp <(cat test-${mod}-${COMP}.$ext | sort) <(cat test-${mod}-${VERSION}.$ext | sort) > /dev/null 2>&1 \
            && echo "  Passed" \
            || echo "  Failed"
    done
done
