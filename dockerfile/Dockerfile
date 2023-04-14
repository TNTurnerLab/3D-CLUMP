FROM continuumio/miniconda3:4.12.0

RUN apt-get update && apt-get install --no-install-recommends -y ca-certificates apt-transport-https libnss-sss git gcc g++ vim less unzip wget tabix && rm -rf /var/lib/apt/lists*

ADD env.yaml /env.yaml
RUN conda env create -f /env.yaml

CMD ["/bin/bash"]
