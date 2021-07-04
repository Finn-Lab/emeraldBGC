FROM continuumio/miniconda3

WORKDIR /docker

COPY environment.yml .
RUN conda env create -f environment.yml

SHELL ["conda", "run", "-n", "bgc", "/bin/bash", "-c"]

RUN python -c "import emeraldbgc"

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bgc", "emeraldbgc"]
