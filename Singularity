BootStrap: docker
From: bioconductor/bioconductor_docker:RELEASE_3_17-R-4.3.0

%files
    ./prepare_data /app/prepare_data
    ./dependencies /app/dependencies

%post
    set -ex
    export DEBIAN_FRONTEND=noninteractive

    # Workaround for headless Java man page directory issue
    mkdir -p /usr/share/man/man1 /usr/share/man/man2

    apt-get update && apt-get install -y --no-install-recommends \
        wget \
        lsb-release

    # randomly copied from the website
    wget https://downloads.skewed.de/skewed-keyring/skewed-keyring_1.5_all_$(lsb_release -s -c).deb
    dpkg -i skewed-keyring_1.5_all_$(lsb_release -s -c).deb

    echo "deb [signed-by=/usr/share/keyrings/skewed-keyring.gpg] https://downloads.skewed.de/apt $(lsb_release -s -c) main" > \
    /etc/apt/sources.list.d/skewed.list

    apt-get update \
        && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        openjdk-11-jre \
        python3-dev \
        python3-pip \
        python-is-python3 \
        r-cran-ggplot2 \
        r-cran-optparse \
        python3-graph-tool \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

    # Install Python and R dependencies
    pip3 install --no-cache-dir -r /app/dependencies/requirements-workflow.txt
    bash /app/dependencies/r-packages-workflow.sh

    # Install C++ and Java dependencies for module detection
    cd /app/prepare_data/workflow/scripts/module_detection

    # -- ClusterONE --
    wget https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar

    # -- LazyFox --
    git clone https://github.com/TimGarrels/LazyFox lazyfoxdir \
        && cd lazyfoxdir \
        && git reset --hard d08f3c084df19bd2a1726159f181bbe3ad6f5bf4 \
        && mkdir build \
        && cd build \
        && cmake .. \
        && make \
        && mv LazyFox ../../LazyFox \
        && cd ../../ \
        && rm -rf lazyfoxdir \
        && chmod +x LazyFox

%environment
    export PATH="/app/prepare_data/workflow/scripts:$PATH"

%runscript
    tail -f /dev/null