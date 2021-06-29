#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Build the documentation from the MAIN_BRANCH or DEV_BRANCH
# and push it to TARGET_BRANCH.
MAIN_BRANCH="main"
TARGET_BRANCH="gh-pages"

mkdir out

# Build the Sphinx documentation
cd sphinx_docs
make html
cd ../

if [ "$GITHUB_BRANCH" = "$MAIN_BRANCH" ]; then
    mv sphinx_docs/build/html/* out/
fi
