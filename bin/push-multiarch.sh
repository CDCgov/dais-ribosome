#!/usr/bin/env bash
# push-multiarch.sh — Build and push linux/amd64 + linux/arm64 manifest to a registry.
#
# Usage:
#   push-multiarch.sh <tag> <repo> <ca_cert_path>
#
#   tag           Image tag, e.g. "1.2.3" or "dev-latest"  (required)
#   repo          Image repository, e.g. "sammysheep/dais-ribosome"  (required)
#   ca_cert_path  Path to a PEM CA certificate to inject at build time via
#                   --secret id=gitlab_ca (required)
#
# Examples:
#   push-multiarch.sh 1.2.3 sammysheep/dais-ribosome ../my-ca.pem
#   push-multiarch.sh dev-latest myorg/dais-ribosome ../my-ca.pem

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TAG="${1:?Usage: push-multiarch.sh <tag> <repo> <ca_cert_path>}"
REPO="${2:?repo required as arg 2 (e.g. sammysheep/dais-ribosome)}"
CA_CERT="${3:?CA cert path required as arg 3 (e.g. ../my-ca.pem)}"

if [[ ! -f "${CA_CERT}" ]]; then
    echo "ERROR: CA cert not found: ${CA_CERT}" >&2
    exit 1
fi

FULL_IMAGE="${REPO}:${TAG}"
PLATFORMS="linux/amd64,linux/arm64"
BUILDER="multiarch-ribosome"

echo "==> Image  : ${FULL_IMAGE}"
echo "==> Platforms: ${PLATFORMS}"

# Ensure buildx builder with multi-platform support exists
if ! docker buildx inspect "${BUILDER}" &> /dev/null; then
    echo "==> Creating buildx builder: ${BUILDER}"
    docker buildx create --name "${BUILDER}" --driver docker-container --bootstrap
fi

docker buildx use "${BUILDER}"

# Assemble build arguments
BUILD_ARGS=(
    buildx build
    --platform "${PLATFORMS}"
    --tag "${FULL_IMAGE}"
    --push
)

BUILD_ARGS+=(--secret "id=gitlab_ca,src=${CA_CERT}")

BUILD_ARGS+=("${REPO_ROOT}")

echo "==> Building and pushing multi-arch manifest..."
docker "${BUILD_ARGS[@]}"

echo "==> Done. Manifest pushed: ${FULL_IMAGE}"
echo "    Platforms: ${PLATFORMS}"
