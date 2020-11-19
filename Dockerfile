FROM python:3.7-slim AS compile-image
RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential gcc
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install MDAnalysis
COPY . /opt/tame
RUN pip install /opt/tame

FROM python:3.7-slim AS build-image
RUN apt-get update -y
RUN apt-get install -y --no-install-recommends libgomp1
COPY --from=compile-image /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
ENTRYPOINT ["tame"]
