

set -ex



samtools --help
samtools view 'https://example.com' 2>&1 | grep 'truncated file.' -q
exit 0
