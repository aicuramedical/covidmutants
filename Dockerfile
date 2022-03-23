FROM python:3.9-slim-buster
ENV PATH=$PATH:/usr/share/miniconda3/bin

VOLUME [ "/output" ]

RUN apt -y update \
   && DEBIAN_FRONTEND=noninteractive apt -y install --no-install-recommends build-essential zlib1g-dev libboost-all-dev libtool git wget emboss \
   && apt clean \
   && rm -rf /var/lib/apt/lists/*

RUN pip3 install "biopython>=1.79" "pandas>=1.4.1" "numpy>=1.18.5" "joblib>=0.14" scipy matplotlib sphinx tabulate

COPY *.py /usr/local/bin/

COPY test/* /output/

RUN git clone https://github.com/TGAC/KAT.git \
  && cd KAT \
  && sed -i 's/\$(top_builddir)\/deps\/boost\/build\/lib/\/usr\/lib\/x86_64-linux-gnu/g' src/Makefile.am \
  && sed -i 's/\$(top_builddir)\/deps\/boost\/build\/lib/\/usr\/lib\/x86_64-linux-gnu/g' lib/Makefile.am \
  && sed -i 's/\$(top_builddir)\/deps\/boost\/build\/lib/\/usr\/lib\/x86_64-linux-gnu/g' tests/Makefile.am \
  && ./autogen.sh \
  && ./configure \
  && make \
  && make install \
  && cd \
  && rm -rf KAT

WORKDIR /output

ENTRYPOINT [ "python3", "/usr/local/bin/oligomutk.py" ]
CMD [ "-h" ]
