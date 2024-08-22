#!/bin/bash

NEW="HEAD"
OLD="${RIBOSOME_OLD_VERSION:-v1.5.4}"

export IFX_LOCAL_PROCS=48

function fix_headers() {
    perl -p -e 's/>.+?\|.+?\|(.+?)$/>$1/' "$1"
}

if [ "$#" == "2" ]; then
    flu_data="cat $1"
    cov_data="cat $2"
    rsv_data="cat $3"
else
    flu_data="fix_headers spec/INFLUENZA.refs"
    cov_data="fix_headers spec/BETACORONAVIRUS.refs"
    rsv_data="fix_headers spec/RSV.refs"
fi

rsv_mod=RSV
flu_mod=INFLUENZA
cov_mod=BETACORONAVIRUS

previous=$(git branch | grep '*' | cut -f2 -d' ')
for tag in "$NEW" "$OLD"; do
    git checkout "$tag" > /dev/null \
        || {
            echo "Git failure could be do to unstashed or uncommited work. Stash/commit and then try again."
            exit 1
        } && ./ribosome install

    for mod in flu cov rsv; do
        this_data="${mod}_data"
        long_mod="${mod}_mod"
        echo "Running '$tag'" \
            && time ./ribosome --module ${!long_mod} <(${!this_data}) \
                "test-${mod}-${tag}.seq" \
                "test-${mod}-${tag}.ins" \
                "test-${mod}-${tag}.del" \
                "test-${mod}-${tag}.gen"
    done

done

for mod in flu cov rsv; do
    for ext in seq ins del gen gen.ins gen.del; do
        echo -n "Checking $mod / $ext" \
            && cmp <(sort "test-${mod}-${OLD}.$ext") <(sort "test-${mod}-${NEW}.$ext") > /dev/null 2>&1 \
            && echo "  Passed" \
            || echo "  Failed"
    done
done

git checkout "$previous"
