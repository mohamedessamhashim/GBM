# Dockerfile
FROM continuumio/miniconda3:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    wget \
    unzip \
    vim \
    libssl-dev \
    libffi-dev \
    python3-dev \
    libxml2-dev \
    libxslt-dev \
    libcurl4-openssl-dev \
    libgsl-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libgfortran5 \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /app

# Copy environment file
COPY environment.yml .

# Create conda environment
RUN conda env create -f environment.yml

# Activate environment
SHELL ["conda", "run", "-n", "gbm_project", "/bin/bash", "-c"]

# Install additional packages that need special handling
RUN conda run -n gbm_project pip install torch-geometric

# Copy project files
COPY . .

# Create data directories
RUN mkdir -p data/raw data/processed data/external results

# Set environment variables
ENV PYTHONPATH=/app/src:$PYTHONPATH
ENV CONDA_DEFAULT_ENV=gbm_project

# Expose Jupyter port
EXPOSE 8888

# Create startup script
RUN echo '#!/bin/bash\n\
conda activate gbm_project\n\
jupyter lab --ip=0.0.0.0 --port=8888 --allow-root --no-browser --NotebookApp.token="" --NotebookApp.password=""\n\
' > /app/start.sh && chmod +x /app/start.sh

# Start Jupyter
CMD ["/app/start.sh"] 