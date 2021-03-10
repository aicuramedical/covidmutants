# covidmutants

# Install requirements
# Download and install docker on MAC

# Mac
https://docs.docker.com/docker-for-mac/install/

# Windows
https://docs.docker.com/docker-for-windows/install/

# Download files
git clone git@github.com:aicuramedical/covidmutants.git
cd covidmutants

# 0. Prepare Environment
Copy fasta files in directry: files

# 1. Build
docker build -t covid .

# 2. Run
# on MAC
docker run -v "$(pwd)"/output:/output -ti covid python3 covid.py -threads 15 -forward TATgCCATTAgTgCAAAgAATAgAgCTCgCAC -reverse GTAATTGGAAcAAGcAAATTcTATGGTGGTTG -amplicon TATGCCATTAGTGCAAAGAATAGAGCTCGCACCGTAGCTGGTGTCTCTATCTGTAGTACTATGACCAATAGACAGTTTCATCAAAAATTATTGAAATCAATAGCCGCCACTAGAGGAGCTACTGTAGTAATTGGAACAAGCAAATTCTATGGTGGTTG -fasta files/gisaid_hcov-19_2021_02_12_21\ B1429.fasta -o /output/output.txt

# 3. Check results in output
ls output

