FROM biopython/biopython
ENV PATH=$PATH:/usr/share/miniconda3/bin

VOLUME [ "/output" ]

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
   && sh Miniconda3-latest-Linux-x86_64.sh -b -p /usr/share/miniconda3 \
   && /usr/share/miniconda3/bin/conda config --add channels defaults \
   && /usr/share/miniconda3/bin/conda config --add channels bioconda \
   && /usr/share/miniconda3/bin/conda config --add channels conda-forge \
   && /usr/share/miniconda3/bin/conda install boost-cpp -y \
   && /usr/share/miniconda3/bin/conda install kat -y \
   && rm -f Miniconda3-latest-Linux-x86_64.sh

COPY covid.py /usr/local/bin
WORKDIR /output

ENTRYPOINT [ "python3", "/usr/local/bin/covid.py" ]
CMD [ "-h" ]