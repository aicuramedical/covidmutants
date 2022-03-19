FROM python:3.9-slim-buster
ENV PATH=$PATH:/usr/share/miniconda3/bin

VOLUME [ "/output" ]

RUN apt -y update \
   && DEBIAN_FRONTEND=noninteractive apt -y install --no-install-recommends wget emboss \
   && apt clean \
   && rm -rf /var/lib/apt/lists/*

RUN pip3 install biopython pandas joblib

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
   && sh Miniconda3-latest-Linux-x86_64.sh -b -p /usr/share/miniconda3 \
   && /usr/share/miniconda3/bin/conda config --add channels bioconda \
   && /usr/share/miniconda3/bin/conda config --add channels conda-forge \
   && /usr/share/miniconda3/bin/conda install kat -y \
   && rm -f Miniconda3-latest-Linux-x86_64.sh

COPY oligomutk.py /usr/local/bin

COPY test/* /output/

WORKDIR /output

ENTRYPOINT [ "python3", "/usr/local/bin/oligomutk.py" ]
CMD [ "-h" ]
