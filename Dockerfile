FROM debian:bookworm-slim AS base

# local apt mirror support
# start every stage with updated apt sources
ARG APT_MIRROR_NAME=
RUN if [ -n "$APT_MIRROR_NAME" ]; then sed -i.bak -E '/security/! s^https?://.+?/(debian|ubuntu)^http://'"$APT_MIRROR_NAME"'/\1^' /etc/apt/sources.list && grep '^deb' /etc/apt/sources.list; fi

RUN apt-get update --allow-releaseinfo-change --fix-missing \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y perl procps locales \
    && sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen \
    && locale-gen \
    && apt clean autoclean \
    && apt autoremove --yes \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/

ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8


FROM base AS builder

# To inject a CA certificate at build time (e.g. for corporate/air-gapped environments):
#   docker build --secret id=gitlab_ca,src=/path/to/cert.pem -t cdcgov/dais-ribosome:test .

COPY . /dais-ribosome

RUN --mount=type=secret,id=gitlab_ca \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y git ca-certificates \
    && if [ -f /run/secrets/gitlab_ca ]; then \
        cp /run/secrets/gitlab_ca /usr/local/share/ca-certificates/gitlab-ca.crt \
        && update-ca-certificates; \
    fi \
    && /dais-ribosome/ribosome install \
    && libp=/dais-ribosome/lib \
    && for i in $(ls "$libp/convert/"|grep -vP 'sam2fasta.pl|delim2fasta.pl|fa2delim.pl|nt2aa.pl');do rm "$libp/convert/$i";done \
    && for i in $(ls "$libp/editMSA/"|grep -vP 'reviseTaxa.pl|codonCorrectStats.pl|stripSequences.pl'); do rm "$libp/editMSA/$i";done \
    && for i in $(ls "$libp/sampling/"|grep -vP 'partitionByField.pl'); do rm "$libp/sampling/$i";done

RUN mkdir /app \
    && mv /dais-ribosome/spec \
    /dais-ribosome/refs \
    /dais-ribosome/bin \
    /dais-ribosome/lib \
    /dais-ribosome/ribosome \
    /dais-ribosome/CHANGELOG.md \
    /dais-ribosome/README.md \
    /dais-ribosome/LICENSE \
    /app/

RUN rm /app/lib/sswsort/bin/ssw_* \
    && ln -s /app/bin/third_party/ssw_Linux_$(uname -m) /app/lib/sswsort/bin/ssw_Linux


FROM base AS final

WORKDIR /app
COPY --from=builder /app /app

ENV IFX_WORK_DIR=/tmp

# Recommended mount point for data volume from host
WORKDIR /data

ENV PATH="/app:/app/lib/sswsort:${PATH}"
