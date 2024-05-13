FROM ubuntu:focal AS base

# local apt mirror support
# start every stage with updated apt sources
ARG APT_MIRROR_NAME=
RUN if [ -n "$APT_MIRROR_NAME" ]; then sed -i.bak -E '/security/! s^https?://.+?/(debian|ubuntu)^http://'"$APT_MIRROR_NAME"'/\1^' /etc/apt/sources.list && grep '^deb' /etc/apt/sources.list; fi

RUN apt-get update --allow-releaseinfo-change --fix-missing \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y git wget ca-certificates \
    && apt clean autoclean \
    && apt autoremove --yes \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/


FROM base AS builder

ARG gitlab_ca
ENV gitlab_ca=${gitlab_ca:-https://docs.cdc.gov/assets/files/CDC-G2-44495b0bcb64fa25e38eea0072929d82.pem}
COPY . /dais-ribosome

RUN wget ${gitlab_ca} \
    && git config --global http.sslCAInfo $(pwd)/$(basename ${gitlab_ca}) \
    && /dais-ribosome/ribosome install \
    && libp=/dais-ribosome/lib \
    && for i in $(ls "$libp/convert/"|grep -vP 'sam2fasta.pl|delim2fasta.pl|fa2delim.pl|nt2aa.pl');do rm "$libp/convert/$i";done \
    && for i in $(ls "$libp/lib/editMSA/"|grep -vP 'reviseTaxa.pl|codonCorrectStats.pl|stripSequences.pl'); do rm "$libp/editMSA/$i";done \
    && for i in $(ls "$libp/sampling/"|grep -vP 'partitionByField.pl'); do rm "$libp/sampling/$i";done \
    && rm -rf /dais-ribosome/workdir /dais-ribosome/lib/sswsort/workdir \
    && ln -s /tmp /dais-ribosome/workdir \
    && ln -s /tmp /dais-ribosome/lib/sswsort/workdir


FROM base as final

COPY --from=builder /dais-ribosome /dais-ribosome

# Recommended mount point for data volume from host
WORKDIR /data

# Export IRMA and LABEL to PATH
ENV PATH "/dais-ribosome:${PATH}"
