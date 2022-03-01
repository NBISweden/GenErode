#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

### Usage: start the script from within the script directory

ncbi_dir="../../testdata/gerp/ncbi_datasets_downloads"

cd ${ncbi_dir}

for i in */
do
  cd ${i}
  unzip *.zip
  cd ..
done