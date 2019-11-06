wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh37.nr_deletions.pathogenic.bed.gz
gunzip GRCh37.nr_deletions.pathogenic.bed.gz
bgzip GRCh37.nr_deletions.pathogenic.bed
tabix GRCh37.nr_deletions.pathogenic.bed.gz
