#! /bin/bash
if [ ! -e methane.pet_small ]; then 
    echo "Downloading PET model from zenodo record"
    wget https://zenodo.org/records/10250171/files/methane.pet_small.zip
    unzip methane.pet_small.zip
fi
echo "Running i-PI"
i-pi input.xml &> log.ipi &
sleep 1
echo "Running driver"
i-pi-py_driver -m pet -o methane.pet_small,ch4.xyz -a pet -u &> log.pet
