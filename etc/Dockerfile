# This container creates an environment similar, but not identical, to
# GitHub Actions.
#
# Suggested build command:
#
#     docker image build --pull --no-cache --tag specter:3.13-slim .
#
# Simplest run command, runs pytest --cov in the container:
#
#     docker run --cpus 2 specter:3.13-slim
#
# Run the container in interactive mode with access to persistent storage:
#
#     docker run --cpus 2 --hostname specter-test --interactive --name specter-test \
#         --rm --tty --volume ${HOME}/Documents/Docker/specter-test:/mnt \
#         specter:3.13-slim /usr/bin/bash
#
#
# "slim" is a minimalist Ubuntu-like base GNU/Linux distribution.
FROM python:3.13-slim

# subversion isn't strictly needed at this time, but if this ever
# turns into a broader test container, then desimodel data might
# be needed, for example.
RUN apt update && \
    apt upgrade --yes && \
    apt install --yes --no-install-recommends git libbz2-dev subversion

# Install base packages;
# install package dependencies;
# install packages needed for testing.
RUN pip install --upgrade --no-cache-dir pip setuptools wheel && \
    pip install --no-cache-dir numba astropy scipy ipython fitsio && \
    pip install --no-cache-dir pytest pytest-cov pytest-astropy pytest-mock

# This directory already exists in the base image, which is handy.
ENV HOME=/home/ubuntu

# Create subdirectory for work.
RUN mkdir -p ${HOME}/git/specter

# Equivalent to cd ${HOME}/git/specter. Interactive
# sessions will also start in this directory, which is nice.
WORKDIR ${HOME}/git/specter

# This is a fairly complicated way of performing a git checkout, but it
# replicates most of the steps performed in a GitHub Actions checkout.
# This also embeds a specific branch and HEAD commit in the container itself.
# In principle the git checkout could be updated with a 'git pull'
# when run in interactive mode.
RUN /usr/bin/git config --global --add safe.directory ${HOME}/git/specter && \
    /usr/bin/git init ${HOME}/git/specter && \
    /usr/bin/git remote add origin https://github.com/desihub/specter && \
    /usr/bin/git config --local gc.auto 0 && \
    /usr/bin/git fetch --no-tags --prune --depth=1 origin && \
    /usr/bin/git checkout numpy-2

# This file sets up a convenient and easily-readable display of versions
# when running pytest.
COPY conftest.py ./

# This command is run in non-interactive sessions, when no other command
# is specified.
CMD [ "pytest", "--cov" ]
