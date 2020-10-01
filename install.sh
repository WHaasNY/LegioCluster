#! /bin/bash


echo "renaming LegioCluster-master folder"
mv LegioCluster-master LegioCluster

echo "making folders and subfolders"
mkdir LegioCluster/input
mkdir LegioCluster/output
mkdir LegioCluster/reads
mkdir LegioCluster/reads/Lpn
mkdir LegioCluster/temp
mkdir LegioCluster/VCF_files
mkdir LegioCluster/VCF_files/Lpn
mkdir LegioCluster/VCF_files/Lpn/ATCC_43290
mkdir LegioCluster/VCF_files/Lpn/D-7158
mkdir LegioCluster/VCF_files/Lpn/F-4185
mkdir LegioCluster/VCF_files/Lpn/Lens
mkdir LegioCluster/VCF_files/Lpn/Paris
mkdir LegioCluster/VCF_files/Lpn/ST23
mkdir LegioCluster/VCF_files/Lpn/All_refs
mkdir LegioCluster/VCF_files/Lpn/Dallas_1E
mkdir LegioCluster/VCF_files/Lpn/F4468
mkdir LegioCluster/VCF_files/Lpn/Lorraine
mkdir LegioCluster/VCF_files/Lpn/Philadelphia_1
mkdir LegioCluster/VCF_files/Lpn/ST42
mkdir LegioCluster/VCF_files/Lpn/Corby
mkdir LegioCluster/VCF_files/Lpn/Detroit-1
mkdir LegioCluster/VCF_files/Lpn/HL06041035
mkdir LegioCluster/VCF_files/Lpn/NCTC11286
mkdir LegioCluster/VCF_files/Lpn/Pontiac
mkdir LegioCluster/VCF_files/Lpn/Toronto-2005

echo "generating links"
ln -s Py_code/add_species.py LegioCluster/add_species.py
ln -s Py_code/LegioCluster_main.py LegioCluster/LegioCluster_main.py
ln -s Py_code/LegioCluster.py LegioCluster/LegioCluster.py

echo "done installing LegioCluster"