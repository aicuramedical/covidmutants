FROM python:3

WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt


#RUN git clone https://github.com/Homebrew/brew ~/.linuxbrew/Homebrew
#RUN mkdir ~/.linuxbrew/bin
#RUN ln -s ../Homebrew/bin/brew ~/.linuxbrew/bin
#RUN eval $(~/.linuxbrew/bin/brew shellenv)

#RUN /root/.linuxbrew/bin/brew install brewsci/bio/kat
COPY . .

# CMD [ "python", "./covid.py" ]

