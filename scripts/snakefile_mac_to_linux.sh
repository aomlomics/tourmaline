# snakefile_mac_to_linux.sh
# After making changes to Snakefile_mac, run this script from main directory 
# to change Mac-specific commands to Linux-specific commands.
cat Snakefile_mac | \
sed "s,gzcat \$line,zcat \$line,g" | \
sed "s,md5 -q,md5sum | awk '{{ print \$1 }}',g" \
> Snakefile_linux
