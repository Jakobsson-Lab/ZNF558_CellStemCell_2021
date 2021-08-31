FROM ubuntu
RUN apt-get update -y
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y bioperl
RUN apt install bedtools
RUN apt-get -y install r-base
RUN apt-get -y install git 
RUN git clone https://github.com/4ureliek/TEanalysis.git
RUN cpan Regexp::Common
# ../TEanalysis/TE-analysis_Shuffle_bed.pl -f peaks_score20.bed -q ./repeatmasker/hg38.fa.out -s bed -r ./GRCh38.p13/hg38.chrom.sizes -e ./GRCh38.p13/GRCh38.p13.genome.fa.gaps.bed 
# tail -n +4 ../for_TEanalysis/hg38.fa.out | awk '{print $5, $6, $7, $0}' OFS="\t" | bedtools intersect -a stdin -b hg38.basic.annotation.noNegative.intron.Intergenic.bed -wa | awk '{print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS="\t" > hg38.fa.intron.intergenic.out

# # docker build --tag bioperl_container -f ./src/Dockerfile .
# # docker run --rm -v /Volumes/My\ Passport/z558/16.08.19/to_docker/:/home -w '/home' --name bioperl_instance -t -i bioperl_container
