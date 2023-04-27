FROM python:3.11-slim-bullseye
ENV PATH=$PATH:/usr/share/miniconda3/bin

VOLUME [ "/output" ]

RUN apt -y update \
   && DEBIAN_FRONTEND=noninteractive apt -y install --no-install-recommends build-essential zlib1g-dev libboost-all-dev libtool git wget emboss \
   && apt clean \
   && rm -rf /var/lib/apt/lists/*

RUN pip3 install "biopython>=1.81" "pandas>=2.0.1" "numpy>=1.24.3" "joblib>=1.2.0"

COPY *.py /usr/local/bin/

COPY test/* /output/

RUN git clone https://github.com/TGAC/KAT.git \
  && cd KAT \
  && sed -i 's/\$(top_builddir)\/deps\/boost\/build\/lib/\/usr\/lib\/x86_64-linux-gnu/g' src/Makefile.am \
  && sed -i 's/\$(top_builddir)\/deps\/boost\/build\/lib/\/usr\/lib\/x86_64-linux-gnu/g' lib/Makefile.am \
  && sed -i 's/\$(top_builddir)\/deps\/boost\/build\/lib/\/usr\/lib\/x86_64-linux-gnu/g' tests/Makefile.am \
  && ./autogen.sh \
  && ./configure --disable-pykat \
  && make \
  && make install \
  && cd \
  && rm -rf KAT

WORKDIR /output

ENTRYPOINT [ "python3", "/usr/local/bin/oligomutk.py" ]
CMD [ "-h" ]
