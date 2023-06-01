FROM continuumio/miniconda3:4.12.0

RUN apt-get update && apt-get install --no-install-recommends -y ca-certificates apt-transport-https libnss-sss git gcc g++ vim less unzip wget tabix && rm -rf /var/lib/apt/lists*

ADD env_slim.yml /env_slim.yml
RUN conda env create -f /env_slim.yml
ENV PATH=/opt/conda/envs/snakemake/bin:$PATH 
RUN echo "export PATH=/opt/conda/envs/snakemake/bin:$PATH" >> ~/.bashrc

CMD ["/bin/bash"]
