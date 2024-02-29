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

## Identify orthologs with proteinortho
```sh
# Before parsing, an error needs to be fixed in the Plasmodium knowlesi genome file:
cd ~/malaria/Resources/Genomes
cat Plasmodium_knowlesi.genome | grep -v '^chromosome' > fixed_Pk.genome

cd ../../../tmp/Prediction/
cp fixed_Pk.gtf ~/malaria/1_Predictions 

# Use gffParse.pl from earlier to print the protein sequences in fasta format for each Plasmodium genome
cd ~/malaria/Resources/Programs/
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_berghei.genome -g  ~/malaria/1_Predictions/genemark.Pb.gtf -b Pb -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_cynomolgi.genome -g  ~/malaria/1_Predictions/genemark.Pc.gtf -b Pc -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/Plasmodium_faciparum.genome -g  ~/malaria/1_Predictions/genemark.Pf.gtf -b Pf -c -p
./gffParse.pl -i ~/malaria/Resources/Genomes/fixed_Pk.genome -g  ~/malaria/1_Predictions/fixed_Pk.gtf -b Pk -c -p
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

# There is a discrepancy in the headers of our fasta files, so you will need to remove the substring from each line starting from the first tab character followed by "length"
# If you do not do this step prior to running proteinortho you will encounter an error in your nohup.out file! 
cd ~/malaria/Parsed
cat Pv.faa | sed 's/\tlength.*//' > new_Pv.faa
cat Pc.faa | sed 's/\tlength.*//' > new_Pc.faa
cat Pk.faa | sed 's/\tlength.*//' > new_Pk.faa
cat Py.faa | sed 's/\tlength.*//' > new_Py.faa
cat Pf.faa | sed 's/\tlength.*//' > new_Pf.faa
cat Tg.faa | sed 's/\tlength.*//' > new_Tg.faa
cat Ht2.faa | sed 's/\tlength.*//' > new_Ht2.faa
cat Pb.faa | sed 's/\tlength.*//' > new_Pb.faa 

mkdir ../Ortho
cp new* ../Ortho

# Run proteinortho like this:
cd ../Ortho
nohup proteinortho6.pl {new_Ht2,new_Pb,new_Pc,new_Pf,new_Pk,new_Pv,new_Py,new_Tg}.faa &

```

# Identifying one-to-one orthologs
```sh
# BUSCO analysis
# Install BUSCO by creating a new conda environment
# More details about proteinortho and installation can be found at: https://busco.ezlab.org/
conda create -n busco -c bioconda busco=5.6.1
conda activate busco

cd
mkdir Busco
cd Busco

# use '-l apicomplexa' to  include only sequences that belong to organisms classified within the taxonomic group Apicomplexa
busco -i ../Ortho/new_Pb.faa -o Pb -m prot -l apicomplexa
busco -i ../Ortho/new_Pc.faa -o Pc -m prot -l apicomplexa
busco -i ../Ortho/new_Pf.faa -o Pf -m prot -l apicomplexa
busco -i ../Ortho/new_Pk.faa -o Pk -m prot -l apicomplexa
busco -i ../Ortho/new_Pv.faa -o Pv -m prot -l apicomplexa
busco -i ../Ortho/new_Py.faa -o Py -m prot -l apicomplexa
busco -i ../Ortho/new_Tg.faa -o Tg -m prot -l apicomplexa
busco -i ../Ortho/new_Ht2.faa -o Ht2 -m prot -l apicomplexa

# move busco results to a new folder 
cd ~/malaria
mkdir Full_Tables
cd Busco
mv Ht2/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Ht2.tsv
mv Pb/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Pb.tsv
mv Pc/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Pc.tsv
mv Pv/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Pv.tsv
mv Py/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Py.tsv
mv Pf/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Pf.tsv
mv Pk/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Pk.tsv
mv Tg/run_apicomplexa_odb10/full_table.csv Full_Tables/full_table_Tg.tsv 

# remove missing or fragmented busco genes 
# some busco genes are duplicated, in which case we want to keep only one sequence ID (sort -k1,1 -u) 
cd ~/malaria
mkdir Complete_Only_Tables 
cd Full_Tables 
cat full_table_Tg.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Tg.tsv
cat full_table_Pf.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Pf.tsv
cat full_table_Py.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Py.tsv
cat full_table_Pk.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Pk.tsv
cat full_table_Pv.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Pv.tsv
cat full_table_Pc.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Pc.tsv
cat full_table_Ht2.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Ht2.tsv
cat full_table_Pb.tsv | awk 'NR>3 {print $1, $2, $3}' | grep -v 'Missing' | grep -v 'Fragmented' | sort -k1,1 -u >> ../Complete_Only_Tables/complete_table_Pb.tsv

# append species name to the sequence ID
cd Complete_Only_Tables
chmod +x add_species.sh
./add_species.sh

# use buscoSearch.py to extract busco genes found in all 8 species 
cd Complete_Only_Tables
python ../Resources/Programs/buscoSearch.py 
cat common_genes_with_species.txt | wc -l #should get 78 genes

# use createFasta.py to generate fasta files for each of the common busco genes 
# this program makes a new directory calls Fasta_files to store all of the generated fasta files 
chmod +x createFasta.py
python createFasta.py 

```
## Alignment 
```sh
# Create a new environment to install clustalo and raxml 
conda create -n clust_and_rax -c bioconda clustalo=1.2.3 raxml=8.2.12
conda activate clust_and_rax

```
