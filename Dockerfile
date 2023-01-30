FROM python:3.10

WORKDIR /api

RUN pip install poetry

COPY ./poetry.lock /api/poetry.lock
COPY ./pyproject.toml /api/pyproject.toml
COPY ./*.py /api/
RUN poetry config virtualenvs.create false
RUN poetry install --no-dev

RUN rm poetry.lock pyproject.toml
RUN pip uninstall --yes poetry

# Install bash (not included by default in alpine image)

# This line is required to be able to use apk
RUN apt-get update
RUN apt-get install bash

COPY ./*.sh /api/
COPY ./config.json /api/

CMD ["bash", "docker_start.sh"]
