Dataset were organized as:

REFID/
    Data/
        FASTQ reads
    manifest.csv

REFID was defined as v[HypoMap Version Number]-[Alphabetical Dataset No.] e.g. v1-07

Scripts:
1. bcl2_formatter.py = convert filenames to bcl2 format for cellranger.
2. NeMO_download = download read files from NeMO manifest
3. SRA_prefetch = prefetch reads from SRA archive
4. SRA_convert = convert .sra files to fastq files