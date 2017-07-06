FROM ibmcom/swift-ubuntu:latest AS build-env

MAINTAINER mtakagi

ADD . /build
WORKDIR /build
RUN swiftc ao.swift -Ounchecked

FROM ibmcom/swift-ubuntu-runtime:latest

COPY --from=build-env /build/ao /usr/local/bin/ao
CMD [ "bash" ]