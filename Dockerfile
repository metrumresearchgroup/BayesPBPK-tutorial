ARG DEFAULT_USER=bayes
FROM alpine:3.15 as rootfs-stage
ARG DEFAULT_USER

# environment
ENV REL=jammy
ENV ARCH=amd64

# install packages
RUN \
  apk add --no-cache \
    bash \
    curl \
    tar \
    tzdata \
    xz \
    git \
    openssh


# RUN git clone git@github.com:timknab/BayesPBPKTutorial.git
RUN git clone https://github.com/metrumresearchgroup/BayesPBPK-tutorial.git BayesPBPKTutorial

# Get base ubuntu files
RUN git clone --branch jammy https://github.com/linuxserver/docker-baseimage-ubuntu.git
# Get vscode files
RUN git clone https://github.com/linuxserver/docker-code-server.git
# Get rocker/rstudio files
RUN git clone https://github.com/rocker-org/rocker-versioned2.git
## Set user
RUN sed -i "s/abc/${DEFAULT_USER}/g" docker-baseimage-ubuntu/root/etc/s6-overlay/s6-rc.d/init-adduser/run
RUN sed -i "s/abc/${DEFAULT_USER}/g" docker-code-server/root/etc/s6-overlay/s6-rc.d/init-code-server/run
RUN sed -i "s/abc/${DEFAULT_USER}/g" docker-code-server/root/usr/local/bin/install-extension
RUN sed -i "s/abc/${DEFAULT_USER}/g" docker-code-server/root/etc/s6-overlay/s6-rc.d/svc-code-server/run


# Change user home
RUN sed -i "s|cp /root/.bashrc /config/.bashrc| cp /root/.bashrc /home/${DEFAULT_USER}/.bashrc|g" docker-code-server/root/etc/s6-overlay/s6-rc.d/init-code-server/run



# grab base tarball
RUN \
  mkdir /root-out && \
  curl -o \
    /rootfs.tar.gz -L \
    https://partner-images.canonical.com/core/${REL}/current/ubuntu-${REL}-core-cloudimg-${ARCH}-root.tar.gz && \
  tar xf \
    /rootfs.tar.gz -C \
    /root-out

# set version for s6 overlay
ARG S6_OVERLAY_VERSION="3.1.2.1"
ARG S6_OVERLAY_ARCH="x86_64"

# add s6 overlay
ADD https://github.com/just-containers/s6-overlay/releases/download/v${S6_OVERLAY_VERSION}/s6-overlay-noarch.tar.xz /tmp
RUN tar -C /root-out -Jxpf /tmp/s6-overlay-noarch.tar.xz
ADD https://github.com/just-containers/s6-overlay/releases/download/v${S6_OVERLAY_VERSION}/s6-overlay-${S6_OVERLAY_ARCH}.tar.xz /tmp
RUN tar -C /root-out -Jxpf /tmp/s6-overlay-${S6_OVERLAY_ARCH}.tar.xz

# add s6 optional symlinks
ADD https://github.com/just-containers/s6-overlay/releases/download/v${S6_OVERLAY_VERSION}/s6-overlay-symlinks-noarch.tar.xz /tmp
RUN tar -C /root-out -Jxpf /tmp/s6-overlay-symlinks-noarch.tar.xz
ADD https://github.com/just-containers/s6-overlay/releases/download/v${S6_OVERLAY_VERSION}/s6-overlay-symlinks-arch.tar.xz /tmp
RUN tar -C /root-out -Jxpf /tmp/s6-overlay-symlinks-arch.tar.xz

# Runtime stage
FROM scratch as ubuntu_base
ARG DEFAULT_USER
COPY --from=rootfs-stage /root-out/ /
ARG BUILD_DATE
ARG VERSION
ARG MODS_VERSION="v3"
LABEL build_version="Linuxserver.io version:- ${VERSION} Build-date:- ${BUILD_DATE}"
LABEL maintainer="TheLamer"

ADD "https://raw.githubusercontent.com/linuxserver/docker-mods/mod-scripts/docker-mods.${MODS_VERSION}" "/docker-mods"

# set environment variables
ARG DEBIAN_FRONTEND="noninteractive"
ENV HOME="/root" \
LANGUAGE="en_US.UTF-8" \
LANG="en_US.UTF-8" \
TERM="xterm" \
S6_CMD_WAIT_FOR_SERVICES_MAXTIME="0" \
S6_VERBOSITY=1 \
S6_STAGE2_HOOK=/docker-mods

# copy sources
COPY --from=rootfs-stage docker-baseimage-ubuntu/sources.list /etc/apt/




RUN \
  echo "**** Ripped from Ubuntu Docker Logic ****" && \
  set -xe && \
  echo '#!/bin/sh' \
    > /usr/sbin/policy-rc.d && \
  echo 'exit 101' \
    >> /usr/sbin/policy-rc.d && \
  chmod +x \
    /usr/sbin/policy-rc.d && \
  dpkg-divert --local --rename --add /sbin/initctl && \
  cp -a \
    /usr/sbin/policy-rc.d \
    /sbin/initctl && \
  sed -i \
    's/^exit.*/exit 0/' \
    /sbin/initctl && \
  echo 'force-unsafe-io' \
    > /etc/dpkg/dpkg.cfg.d/docker-apt-speedup && \
  echo 'DPkg::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' \
    > /etc/apt/apt.conf.d/docker-clean && \
  echo 'APT::Update::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' \
    >> /etc/apt/apt.conf.d/docker-clean && \
  echo 'Dir::Cache::pkgcache ""; Dir::Cache::srcpkgcache "";' \
    >> /etc/apt/apt.conf.d/docker-clean && \
  echo 'Acquire::Languages "none";' \
    > /etc/apt/apt.conf.d/docker-no-languages && \
  echo 'Acquire::GzipIndexes "true"; Acquire::CompressionTypes::Order:: "gz";' \
    > /etc/apt/apt.conf.d/docker-gzip-indexes && \
  echo 'Apt::AutoRemove::SuggestsImportant "false";' \
    > /etc/apt/apt.conf.d/docker-autoremove-suggests && \
  mkdir -p /run/systemd && \
  echo 'docker' \
    > /run/systemd/container && \
  echo "**** install apt-utils and locales ****" && \
  apt-get update && \
  apt-get install -y \
    apt-utils \
    locales && \
  echo "**** install packages ****" && \
  apt-get install -y \
    curl \
    gnupg \
    jq \
    netcat \
    tzdata && \
  echo "**** generate locale ****" && \
  locale-gen en_US.UTF-8 && \
  echo "**** create default user and make our folders ****" && \
  useradd -u 911 -U -s /bin/false ${DEFAULT_USER} && \
  usermod -G users ${DEFAULT_USER} && \
  mkdir -p \
    /app \
    /config \
    /home/${DEFAULT_USER} \
    /defaults && \
  chown -R ${DEFAULT_USER}:${DEFAULT_USER} /home/${DEFAULT_USER} && \
  chmod +x /docker-mods && \
  echo "**** cleanup ****" && \
  apt-get autoremove && \
  apt-get clean && \
  rm -rf \
    /tmp/* \
    /var/lib/apt/lists/* \
    /var/tmp/* \
    /var/log/*


COPY --from=rootfs-stage docker-baseimage-ubuntu/root/ /

ENTRYPOINT ["/init"]


FROM ubuntu_base as vscode_base
# COPY OVER PROJECT FILES
COPY --from=rootfs-stage BayesPBPKTutorial /home/${DEFAULT_USER}/BayesPBPK
ENV DEFAULT_WORKSPACE=/home/${DEFAULT_USER}/BayesPBPK
# set version label
ARG BUILD_DATE
ARG VERSION
ARG CODE_RELEASE
LABEL build_version="Linuxserver.io version:- ${VERSION} Build-date:- ${BUILD_DATE}"
LABEL maintainer="TimKnab"

# environment settings
ARG DEBIAN_FRONTEND="noninteractive"
ENV HOME="/home/${DEFAULT_USER}"

RUN \
  echo "**** install runtime dependencies ****" && \
  apt-get update && \
  apt-get install -y \
    git \
    jq \
    libatomic1 \
    nano \
    net-tools \
    netcat \
    sudo && \
  echo "**** install code-server ****" && \
  if [ -z ${CODE_RELEASE+x} ]; then \
    CODE_RELEASE=$(curl -sX GET https://api.github.com/repos/coder/code-server/releases/latest \
      | awk '/tag_name/{print $4;exit}' FS='[""]' | sed 's|^v||'); \
  fi && \
  mkdir -p /app/code-server && \
  curl -o \
    /tmp/code-server.tar.gz -L \
    "https://github.com/coder/code-server/releases/download/v${CODE_RELEASE}/code-server-${CODE_RELEASE}-linux-amd64.tar.gz" && \
  tar xf /tmp/code-server.tar.gz -C \
    /app/code-server --strip-components=1 && \
  echo "**** clean up ****" && \
  apt-get clean && \
  rm -rf \
    /config/* \
    /tmp/* \
    /var/lib/apt/lists/* \
    /var/tmp/*


# add local files
COPY --from=rootfs-stage docker-code-server/root /

## Install extensions
RUN /app/code-server/bin/code-server --extensions-dir /config/extensions --install-extension julialang.language-julia
RUN /app/code-server/bin/code-server --extensions-dir /config/extensions --install-extension synedra.auto-run-command
RUN /app/code-server/bin/code-server --extensions-dir /config/extensions --install-extension Ikuyadeu.r


# ports and volumes
EXPOSE 8443
RUN cp /init /init_lsio


FROM vscode_base as rbase

ENV R_VERSION=4.2.1
ENV R_HOME=/usr/local/lib/R
ENV TZ=Etc/UTC


COPY --from=rootfs-stage rocker-versioned2/scripts/install_R_source.sh /rocker_scripts/install_R_source.sh
RUN /rocker_scripts/install_R_source.sh



ENV CRAN=https://packagemanager.rstudio.com/cran/__linux__/jammy/latest
ENV LANG=en_US.UTF-8

COPY --from=rootfs-stage rocker-versioned2/scripts /rocker_scripts

RUN /rocker_scripts/setup_R.sh

FROM rbase as rstudio
ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2022.07.2+576
ARG DEFAULT_USER
ENV DEFAULT_USER=bayes
ENV PANDOC_VERSION=default
ENV QUARTO_VERSION=default

RUN /rocker_scripts/install_rstudio.sh && \
    /rocker_scripts/install_pandoc.sh && \
    /rocker_scripts/install_quarto.sh && \
    /rocker_scripts/install_julia.sh && \
    cp /init_lsio /init && \
    echo "auth-minimum-user-id=911" >> /etc/rstudio/rserver.conf && \
    rm -rf /home/${DEFAULT_USER}/R

EXPOSE 8787



## Install R Packages

RUN rm -rf /home/${DEFAULT_USER}/R && \
    chown -R ${DEFAULT_USER}:${DEFAULT_USER} /home/${DEFAULT_USER} && \
    R -e "install.packages('rstudioapi')" && \
    R -e "install.packages('languageserver')"
COPY .Rprofile /home/${DEFAULT_USER}/.Rprofile
RUN chown ${DEFAULT_USER}:${DEFAULT_USER} /home/${DEFAULT_USER}/.Rprofile 

FROM rstudio as final
COPY settings.json /config/data/Machine/settings.json
COPY settings_user.json /config/data/User/settings.json

WORKDIR /home/${DEFAULT_USER}/BayesPBPK/script
RUN git clone https://github.com/metrumresearchgroup/Torsten && \
    julia --project=/home/${DEFAULT_USER}/BayesPBPK -e 'import Pkg; Pkg.instantiate()'  && \
    chown -R ${DEFAULT_USER}:${DEFAULT_USER} /home/${DEFAULT_USER} && \
    rm -rf /home/${DEFAULT_USER}/BayesPBPK/.Rproj.user && \
    rm -rf /init_lsio && \
    rm -rf /rocker_scripts && \
    rm -rf /home/${DEFAULT_USER}/BayesPBPK/script/Torsten/.git && \
    rm -rf /home/${DEFAULT_USER}/BayesPBPK/script/Torsten/example-models && \
    rm -rf /home/${DEFAULT_USER}/BayesPBPK/script/Torsten/.gitignore && \
    rm -rf /home/${DEFAULT_USER}/BayesPBPK/script/Torsten/.github

ENV PASSWORD=metrumrg
RUN sed -i 's|AUTH="password"|AUTH="none"|g'  /etc/s6-overlay/s6-rc.d/svc-code-server/run
RUN echo "script/Torsten/" >> /home/${DEFAULT_USER}/BayesPBPK/.gitignore
ENV SUDO_PASSWORD='metrumrg'

## Install additional packages
RUN apt update && \
    apt install -y \
      zlib1g-dev \
      libxml2-dev \
      libfontconfig1-dev \
      ncdu \
      pkg-config \
      curl

WORKDIR /home/${DEFAULT_USER}/BayesPBPK
RUN R -e 'install.packages("renv"); install.packages("jsonlite"); renv::activate("/home/'${DEFAULT_USER}'/BayesPBPK"); renv::install("jsonlite"); renv::restore(); renv::install("languageserver")'
RUN mkdir /data_store 
RUN chown -R ${DEFAULT_USER}:${DEFAULT_USER} /data_store
# Run with  docker run -it --rm  -p 8443:8443 -p 8787:8787 bayespbpk:latest


