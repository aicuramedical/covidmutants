# covidmutants

# Install requirements
# Download and install docker on MAC

# Mac
https://docs.docker.com/docker-for-mac/install/

# Windows
https://docs.docker.com/docker-for-windows/install/

# 1. Build
docker build -t covid .

# 2. Run
docker run -ti covid python covid.py -threads 15 -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG -fasta files/gisaid_hcov-19_2021_02_12_21\ B1429.fasta

