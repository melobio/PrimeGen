FROM 10.49.2.32:8050/python:3.12-slim-bookworm
LABEL maintainer="chenjialang@genomics.cn"

ENV LANG=C.UTF-8
ENV PATH="$PATH:/app/bin"

COPY ./depends/sources.list /etc/apt/sources.list
COPY ./depends/pip.conf /root/.pip/
COPY ./depends/datasets /app/bin/
COPY ./depends/taxdump.tar.gz /tmp/
COPY ./depends/taxonkit /app/bin/
COPY ./depends/bedtools/* /app/bin/
COPY ./depends/blast/* /app/bin/

RUN chmod +x /app/bin/datasets \
    && chmod +x /app/bin/taxonkit \
    && mkdir /root/.taxonkit \
    && tar -zxvf /tmp/taxdump.tar.gz -C /root/.taxonkit 
    

RUN  apt-get clean \
    && apt-get update \
    && apt-get install -f -y --no-install-recommends \
    ca-certificates \
    curl \
    dirmngr \
    gnupg \
    xz-utils \
    vim \
    gcc \
    build-essential \
    libsasl2-dev python3-dev libldap2-dev libssl-dev \
    && apt-get autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*

# 容器内工作目录
WORKDIR /app/src

# 第三方依赖包
COPY ./lib /app/lib
RUN pip3 install /app/lib/*.whl
COPY ./requirements.txt /app/requirements.txt
RUN pip3 install -r /app/requirements.txt

ENV HOST="0.0.0.0"
ENV PORT=8081
ENV WEB_CONCURRENCY=2

EXPOSE $PORT

# 项目代码
COPY ./src /app/src

VOLUME [ "/reference_data", "/app/logs" ]

# COPY ./datasets /bin/datasets

CMD uvicorn main:app --host $HOST --port ${PORT}
