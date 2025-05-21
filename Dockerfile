FROM python:3.9

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

RUN pip install --upgrade pip
RUN pip install pysam
RUN pip install pandas
RUN pip install scipy

COPY *.py /scripts/
COPY pipseq_barcode_sequence.txt /scripts/