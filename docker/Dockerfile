FROM continuumio/miniconda3

WORKDIR /docker

COPY environment.yml .
COPY get_ips_slim.sh .
RUN conda env create -f environment.yml
RUN conda run --no-capture-output -n bgc bash get_ips_slim.sh

SHELL ["conda", "run", "-n", "bgc", "/bin/bash", "-c"]

RUN python -c "import emeraldbgc"

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bgc", "emeraldbgc"]
