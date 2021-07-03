FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "bgc", "/bin/bash", "-c"]

# Demonstrate the environment is activated:
RUN python -c "import emeraldbgc"

# The code to run when container is started:
#COPY run.py .
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bgc", "emeraldbgc"]
