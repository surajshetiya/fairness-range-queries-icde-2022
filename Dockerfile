# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

# Update apps on the base image
RUN apt-get -y update && apt-get install -y

# Set up timezone for boost
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata

# RUN timedatectl set-timezone America/Chicago

RUN echo "America/Chicago" | tee /etc/timezone
RUN ln -s -f /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata

# Install the Clang compiler
RUN apt-get -y install clang libboost-all-dev cmake lldb