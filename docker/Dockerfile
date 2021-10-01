FROM continuumio/miniconda3

LABEL base_image="continuumio/miniconda3"
LABEL version="1"
LABEL about.summary="emeraldBGC - SMBGC detection tool -"
LABEL about.license="SPDX:Apache-2.0"
LABEL about.tags="SMBGC BGC"
LABEL about.home="https://github.com/Finn-Lab/emeraldBGC"
LABEL software="emeraldbgc"

WORKDIR /docker

COPY environment.yml .
COPY get_ips_slim.sh .

RUN conda env create -f environment.yml
RUN conda run --no-capture-output -n bgc bash get_ips_slim.sh

SHELL ["conda", "run", "-n", "bgc", "/bin/bash", "-c"]

RUN python -c "import emeraldbgc"

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bgc", "emeraldbgc"]
