FROM biopython/biopython

WORKDIR /usr/src/app

VOLUME [ "/output" ]
# COPY requirements.txt ./
# RUN pip install --no-cache-dir -r requirements.txt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN sh Miniconda3-latest-Linux-x86_64.sh -b
RUN /root/miniconda3/bin/conda config --add channels defaults
RUN /root/miniconda3/bin/conda config --add channels bioconda
RUN /root/miniconda3/bin/conda config --add channels conda-forge

RUN /root/miniconda3/bin/conda install bwa -y
RUN /root/miniconda3/bin/conda install kat -y
#RUN git clone https://github.com/Homebrew/brew ~/.linuxbrew/Homebrew
#RUN mkdir ~/.linuxbrew/bin
#RUN ln -s ../Homebrew/bin/brew ~/.linuxbrew/bin
#RUN eval $(~/.linuxbrew/bin/brew shellenv)

#RUN /root/.linuxbrew/bin/brew install brewsci/bio/kat
COPY . .

# CMD [ "python", "./covid.py" ]

