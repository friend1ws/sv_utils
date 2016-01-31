#! /usr/bin/env bash

# create corresponding table for GRCh name and UCSC name
rm -rf GCF_000001405.13.assembly.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.13.assembly.txt
python make_ucsc_grch.py GCF_000001405.13.assembly.txt > grch2ucsc.txt

# for GRCh38 (hg38)
# ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.31.assembly.txt

rm -rf refGene.txt.gz
rm -rf ensGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz
 
echo "python listGene.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refGene.bed"
python listGene.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refGene.bed

echo "python listGene.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensGene.bed"
python listGene.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensGene.bed

echo "python listExon.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refExon.bed"
python listExon.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refExon.bed

echo "python listExon.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensExon.bed"
python listExon.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensExon.bed

echo "python listJunc.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refJunc.bed"
python listJunc.py refGene.txt.gz ref | sort -k1,1 -k2,2n -k3,3n - > refJunc.bed

echo "python listJunc.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensJunc.bed"
python listJunc.py ensGene.txt.gz ens | sort -k1,1 -k2,2n -k3,3n - > ensJunc.bed

echo "python listCoding.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n > refCoding.bed"
python listCoding.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n > refCoding.bed

# compressing  and indexing
echo "bgzip -f refGene.bed"
bgzip -f refGene.bed

echo "bgzip -f ensGene.bed"
bgzip -f ensGene.bed

echo "bgzip -f refExon.bed"
bgzip -f refExon.bed

echo "bgzip -f ensExon.bed"
bgzip -f ensExon.bed

echo "bgzip -f refJunc.bed"
bgzip -f refJunc.bed

echo "bgzip -f ensJunc.bed"
bgzip -f ensJunc.bed

echo "bgzip -f refCoding.bed"
bgzip -f refCoding.bed

echo "tabix -f -p bed refGene.bed.gz"
tabix -f -p bed refGene.bed.gz 

echo "tabix -f -p bed ensGene.bed.gz"
tabix -f -p bed ensGene.bed.gz 

echo "tabix -f -p bed refExon.bed.gz"
tabix -f -p bed refExon.bed.gz 

echo "tabix -f -p bed ensExon.bed.gz"
tabix -f -p bed ensExon.bed.gz 

echo "tabix -f -p bed refJunc.bed.gz"
tabix -f -p bed refJunc.bed.gz

echo "tabix -f -p bed ensJunc.bed.gz"
tabix -f -p bed ensJunc.bed.gz

echo "tabix -f -p bed refCoding.bed.gz"
tabix -f -p bed refCoding.bed.gz

# simple repeat
rm -rf simpleRepeat.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
zcat simpleRepeat.txt.gz | cut -f 2- > simpleRepeat.bed

echo "bgzip -f simpleRepeat.bed"
bgzip -f simpleRepeat.bed

echo "tabix -f -p bed simpleRepeat.bed.gz"
tabix -f -p bed simpleRepeat.bed.gz


