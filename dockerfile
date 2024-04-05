# Use the latest Ubuntu LTS as the base image
FROM ubuntu:22.04

# Update the package lists and install required dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    build-essential \
    gcc \
    g++ \
    python3.12 \
    python3.12-dev \
    python3-pip \
    cuda-toolkit-11-6 \
    libcudnn8 \
    && rm -rf /var/lib/apt/lists/*

# Set the default Python version to 3.12
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 1
RUN update-alternatives --config python3

# Install the latest version of pip
RUN python3 -m pip install --upgrade pip

# Set the working directory
WORKDIR /app

# Copy your application code into the container
COPY . .

# Install Python dependencies
RUN pip3 install -r requirements.txt

# Set the default command to run your application
CMD ["python3", "main.py"]