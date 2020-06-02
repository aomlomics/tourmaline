#!/bin/bash

# Usage: scripts/snakefile_mac_to_linux.sh
#
# From the main directory run this script after making changes to Snakefile_mac
# to change Mac-specific commands to Linux-specific commands in Snakefile_linux.

cat Snakefile_mac | \
sed "s,gzcat \$line,zcat \$line,g" | \
sed "s,md5 -q,md5sum | awk '{{ print \$1 }}',g" \
> Snakefile_linux
