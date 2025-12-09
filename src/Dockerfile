FROM analysiscommon/genesis_wdl:v1.5

MAINTAINER Seung Hoan Choi <seuchoi@bu.edu>

USER root

RUN sudo wget https://bit.ly/3rmveIc -O rpackages4_1_3_aou.tar.gz
RUN tar -xvf rpackages4_1_3_aou.tar.gz

RUN sudo wget https://bioconductor.statistik.tu-dortmund.de/packages/3.11/bioc/src/contrib/GENESIS_2.18.0.tar.gz
RUN sudo R CMD INSTALL GENESIS_2.18.0.tar.gz

RUN sudo git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git
RUN sudo git clone https://github.com/broadinstitute/phewas_genetic_analysis.git
RUN sudo git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
RUN sudo git clone --branch AOU_rareV_code https://github.com/seuchoi/meta_transcript.git

## how to build the docker
#docker rmi -f e0137a0cac4b
#docker build --no-cache -t us.gcr.io/aou-project-385720/meta_transcript_2.18.0:10282025 ./
#docker push us.gcr.io/aou-project-385720/meta_transcript_2.18.0:10282025
#docker tag e0137a0cac4b genesis_meta_2.18.0:Sep172024
#docker tag genesis_meta_2.18.0:Sep172024 us.gcr.io/aou-project-385720/genesis_meta_2.18.0:Sep172024
#docker push us.gcr.io/aou-project-385720/genesis_meta_2.18.0:Sep172024
#docker tag us.gcr.io/aou-project-385720/genesis_2.18.0:feb132024 schoi/genesis_2.18.0:feb132024
#docker push --disable-content-trust schoi/genesis_2.18.0:feb132024
