FROM ubuntu:20.04 as builder

# Locale
ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# ensure any pipes fail correctly, may impact other things
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# work in an area that will get discarded
WORKDIR /tmp/build

# apt stuff
# hadolint ignore=DL3008
RUN apt-get -yq update \
&& apt-get -yq install --no-install-recommends software-properties-common \
&& add-apt-repository ppa:deadsnakes/ppa \
&& apt-get -yq install --no-install-recommends python3.9 \
&& update-alternatives --install /usr/bin/python python /usr/bin/python3.9 1

# not in final image
# hadolint ignore=DL3008
RUN apt-get -yq update \
&& apt-get -yq install --no-install-recommends \
    build-essential \
    python3.9-dev \
    python3.9-distutils \
    python3.9-venv \
    curl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libmagic-dev \
    default-jdk \
&& curl -sSL --retry 10 https://bootstrap.pypa.io/get-pip.py > get-pip.py \
&& python3.9 get-pip.py

# base environment
ENV OPT /opt/wsi-t78
ENV VIRTUAL_ENV=$OPT/venv
RUN python3.9 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"


#Install smart-phase
RUN curl -sSL --retry 10 https://github.com/paulhager/smart-phase/archive/refs/tags/v1.2.0.tar.gz > sphase.tar.gz \
    && tar -zxf sphase.tar.gz
WORKDIR /tmp/build//smart-phase-1.2.0/
RUN bash compile.sh \
    && cp smartPhase.jar $OPT/

WORKDIR /tmp/build

COPY ./python .

# deploy properly
# hadolint ignore=DL3013
RUN pip install --no-cache-dir ./


########################## FINAL IMAGE ##########################
FROM ubuntu:20.04

# Locale
ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# apt stuff
# hadolint ignore=DL3008
RUN apt-get -yq update \
&& apt-get -yq install --no-install-recommends software-properties-common libmagic-dev zlib1g curl \
&& add-apt-repository ppa:deadsnakes/ppa \
&& apt-get -yq install --no-install-recommends python3.9 default-jre\
# make sure all security patches are applied
&& apt-get -yq update && apt-get install -yq --no-install-recommends unattended-upgrades \
&& unattended-upgrade -d -v \
&& apt-get remove -yq unattended-upgrades \
&& apt-get autoremove -yq \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/* \
&& update-alternatives --install /usr/bin/python python /usr/bin/python3.9 1

ENV OPT /opt/wsi-t78
ENV VIRTUAL_ENV=$OPT/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

COPY --from=builder $OPT $OPT

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu
USER ubuntu

RUN casmsmartphase generate-bed --version

WORKDIR /home/ubuntu

CMD ["java","-jar", "/opt/wsi-t78/smartPhase.jar"]
