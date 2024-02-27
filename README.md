## Getting Started 
```sh
ssh Server
mkdir malaria 
cd malaria 
cp ../../../binp29/Data/plasmodiumGenes.tgz
tar -xvzf plasmodiumGenes.tgz
```

## Gene Prediction of Plasmodium Data
```sh
# check permissions of GeneMark key
# should look like this: 
ls -l .gm_key 

# once you have the correct permissions, run GeneMark for Plasmodium berghei genome
# since berghei has contigs smaller than 50,000 (default) you will need to set --min_config to 10000 
gmes_petap.pl --ES --min_contig 10000 --sequence Plasmodium_berghei.genome

# add results to Predictions folder on server 
cp genemark.gtf ../../../tmp/Predictions/genemark.Pb.gtf

# since it takes a while to run GeneMark for each genome, copy the other predictions from Predictions folder on the server
mkdir 1_Predictions 
mv genemark.gtf 1_Predictions/genemark.Pb.gtf

cd ../../../tmp/Prediction
cp Plasmodium_falciparum/genemark.gtf ~/malaria/genemark.Pf.gtf 
cp genemark.Pk.gtf ~/malaria/1_Predictions/genemark.Pk.gtf
cp genemark.Pc.gtf ~/malaria/1_Predictions/genemark.Pc.gtf
cp genemark.Pv.gtf ~/malaria/1_Predictions/genemark.Pv.gtf
cp genemark.Py.gtf ~/malaria/1_Predictions/genemark.Py.gtf

# download Toxoplasma gff-file from the course website 
cd 
cp ../../../resources/binp29/Data/malaria/Tg.gff.gz ~/malaria/1_Predictions/ 
gunzip Tg.gff.gz
```

## Processing of H. tartakovskyi data
```sh
# download and unzip scaffold file  
cd ~
cp ../../../resources/binp29/Data/malaria/Haemoproteus_tartakovskyi.raw.genome.gz ~/malaria
cd ~/malaria
gunzip Haemoproteus_tartakovskyi.raw.genome.gz

# download removeScaffold.py program from course server
cd malaria
mkdir Programs
cp ../../../../resources/binp29/Data/removeScaffold.py ~/malaria/Programs
cd Programs
chmod +xrw removeScaffold.py

# remove scaffolds above GC-content of 25% to remove bird reads (determined from https://portal.research.lu.se/en/publications/the-genome-of-haemoproteus-tartakovskyi-and-its-relationship-to-h)
# remove scaffolds that are less than 3000 nucleotides 
./removeScaffold.py ~/malaria/Resources/Genomes/Haemoproteus_tartakovskyi.raw.genome 25 ~/malaria/Ht.genome 3000

# check number of scaffolds in clean genome file 
cat Ht.genome | grep '^>' | wc -l # should get 1010

# make a gene prediction from the new genome file
gmes_petap.pl --ES --min_contig 10000 --sequence Ht.genome

```

## Blastp 
```sh
# run blastp to decide if the query sequence orginates from the avian host 
# WARNING: do NOT use more than 10 threads
# I used 40 since I ran this at 1 AM :) 
blastp -query Ht2.faa -db SwissProt -evalue 1e-10 -out Ht_blast_results -num_descriptions 10 -num_alignments 5 -num_threads 40

```
