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
cd ~/malaria
cat Ht.genome | grep '^>' | wc -l # should get 1010

# make a gene prediction from the new genome file
gmes_petap.pl --ES --min_contig 10000 --sequence Ht.genome

```

## Blastp 
```sh
# run blastp to search for genes that are from avian origin
# WARNING: do NOT use more than 10 threads
# I used 40 since I ran this at 1 AM :) 
cd ~/malaria
blastp -query Ht2.faa -db SwissProt -evalue 1e-10 -out Ht_blast_results -num_descriptions 10 -num_alignments 5 -num_threads 40

# Use the taxonomy for the top hit in the BLAST output to decide if the query sequence orginates from the avian host
# the files are approximately 4GB, so instead of copying these files you can create a symlink:
ln -s /resources/binp29/Data/malaria/taxonomy.dat taxonomy.dat
ln -s /resources/binp29/Data/malaria/uniprot_sprot.dat uniprot_sprot.dat
ls -l

# You can also download the taxonomy.dat file from NCBI like this:
# wget ftp://ftp.ebi.ac.uk/pub/databases/taxonomy/taxonomy.dat
# and the SwissProt dat file:
# wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz


# copy the datParser.py script to your Programs folder 
cd ../../../../resources/binp29/Data
cp datParser.py ~/malaria/Resources/Programs/

# Identify scaffolds from genes with avian origin
cd malaria
python ~/malaria/Resources/Programs/datParser.py Ht_blast_results Ht.fna taxonomy.dat uniprot_sprot.dat > scaffolds.txt

# You should get an empty file, which indicates that all avian scaffolds were removed with removeScaffold.py
# I checked to make sure this was correct by using Lovisa's .fna file, which should result in 6 scaffolds:
python ~/malaria/Resources/Programs/datParser.py Ht_blast_results Htartakovskyi.fna taxonomy.dat uniprot_sprot.dat > scaffolds_lovisa.txt

```

## Identify orthologs with proteinortho and BUSCO
```sh
# Use gffParse.pl from earlier to print the protein sequences in fasta format for each Plasmodium genome
cd ~/malaria/Resources/Programs/
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_berghei.genome -g  ~/malaria/1_Predictions/genemark.Pb.gtf -b Pb -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_cynomolgi.genome -g  ~/malaria/1_Predictions/genemark.Pc.gtf -b Pc -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_faciparum.genome -g  ~/malaria/1_Predictions/genemark.Pf.gtf -b Pf -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_knowlesi.genome -g  ~/malaria/1_Predictions/genemark.Pk.gtf -b Pk -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_vivax.genome -g  ~/malaria/1_Predictions/genemark.Pv.gtf -b Pv -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_yoelii.genome -g  ~/malaria/1_Predictions/genemark.Py.gtf -b Py -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Toxoplasma_gondii.genome -g  ~/malaria/1_Predictions/Tg.gff -b Tg -c -p

# Move all of the output to a directory called Parsed for easy access
mkdir ../../Parsed
mv *.log ../../Parsed/
mv *.faa ../../Parsed/
mv *.fna ../../Parsed/

# Install proteinortho with howebrew
# Installation instructions for homebrew can be found at: https://docs.brew.sh/Installation
git clone https://github.com/Homebrew/brew homebrew
eval "$(homebrew/bin/brew shellenv)"
brew update --force --quiet
chmod -R go-w "$(brew --prefix)/share/zsh"

# More details about proteinortho and installation can be found at: https://www.bioinf.uni-leipzig.de/Software/proteinortho/
brew install proteinortho=6.3.1
 
cd ~/Parsed # if you are not already in this directory 

# Run proteinortho like this:
nohup proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &

# BUSCO analysis
# Install BUSCO with homebrew
# More details about proteinortho and installation can be found at: https://busco.ezlab.org/
brew install BUSCO=5.6.1 

cd
mkdir BUSCO 
cd BUSCO

# use '-l apicomplexa' to  include only sequences that belong to organisms classified within the taxonomic group Apicomplexa
busco -i ../Parsed/Pb.faa -o Pb -m prot -l apicomplexa
busco -i ../Parsed/Pc.faa -o Pc -m prot -l apicomplexa
busco -i ../Parsed/Pf.faa -o Pf -m prot -l apicomplexa
busco -i ../Parsed/Pk.faa -o Pk -m prot -l apicomplexa
busco -i ../Parsed/Pv.faa -o Pv -m prot -l apicomplexa
busco -i ../Parsed/Py.faa -o Py -m prot -l apicomplexa
busco -i ../Parsed/Tg.faa -o Tg -m prot -l apicomplexa

```

