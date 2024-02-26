## Getting Started 
ssh Server
mkdir malaria 
cd malaria 
cp ../../../binp29/Data/plasmodiumGenes.tgz
tar -xvzf plasmodiumGenes.tgz

## Gene Prediction of Plasmodium Data
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
cp #Pv
cp genemark.Py.gtf ~/malaria/1_Predictions/genemark.Py.gtf

# download Toxoplasma gff-file from the course website 
# need to do this 

## Processing of H. tartakovskyi data
# download scaffold file  
cd ~
cp ../../../resources/binp29/Data/malaria/Haemoproteus_tartakovskyi.raw.genome.gz ~/malaria
cd ~/malaria
gunzip Haemoproteus_tartakovskyi.raw.genome.gz
