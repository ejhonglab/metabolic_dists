#!/usr/bin/env bash

# NOTE: originally ran most of this from ~/src/metabolike, when the config files were
# still there
# TODO test/update now that i moved config to this repo (metabolic_dists)

#######################################################################################
# https://docs.docker.com/engine/install/ubuntu/
#######################################################################################

sudo apt-get remove docker docker-engine docker.io containerd runc

sudo apt-get update
sudo apt-get install \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

sudo mkdir -m 0755 -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin


#######################################################################################
# https://metabolike.readthedocs.io/en/latest/usage/installation/
#######################################################################################
# (assumes this exists. it is where i conventionally keep all my git repos.)
cd ~/src

git clone --depth 1 https://github.com/y1zhou/metabolike


# my own modifications
# (i used mamba here, as it is generally faster)
conda create -n "metabolike" python=3.10

conda activate metabolike

cd metabolike
pip install -e .

#######################################################################################
# https://metabolike.readthedocs.io/en/latest/usage/getting_started/
#######################################################################################

# (assumes docker-compose.yaml and .env are already in place)
sudo docker compose up -d

# (assumes metabolike-setup.yaml is in place)
# NOTE: I found I needed to add the --no-create-db, and even though I didn't manually
# create a db, it still worked and persisted across `docker compose stop`
metabolike setup -c metabolike-setup.yaml --no-create-db
