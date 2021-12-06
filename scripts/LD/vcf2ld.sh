#!bin/bash
##################################################################
###   Pairwise SNP Linkage Disequilibrium Satistics From MSA   ###
##################################################################

# this script requires bcftools and tomahawk to be installed and symlinked
# it will produce a diploid .bcf and tomahawk intermediate files (.twk, .two) and the final pairwise ld table (.txt)

# usage:
# bash vcf2ld.sh <(path/to/)diploid.vcf> <(output/path/)prefix>

# get path/diploid.vcf as first argument
input=$1

diploid_check=$(tail -2 $input | head -1 | awk '{print $10}' | grep -o "/")
if [[ -z $diploid_check ]] ; then
	echo "looks like you're using a haploid vcf, please diploidize it"
	exit
fi

# get output prefix (and path, if specifying different output directory)
# if ignored, defaults to msa prefix and directory
prefix=$2
if [[ -z $prefix ]] ; then
	prefix=$(echo ${input} | rev | awk 'BEGIN{FS="."}{for (i=2; i<NF; i++) printf $i "."; print $NF}' | rev)
fi

# convert diploid vcf to binary
bcftools view ${input} -O b -o ${prefix}.bcf

# run tomahawk
tomahawk import -f -i ${prefix}.bcf -o ${prefix}_LD
tomahawk calc -u -i ${prefix}_LD.twk -o ${prefix}_LD -t 6 -r 0
tomahawk view -i ${prefix}_LD.two > ${prefix}_LD_table.txt
tail -n +14 ${prefix}_LD_table.txt > ${prefix}_LD_table_clean.txt
head -13 ${prefix}_LD_table.txt > ${prefix}_tomahawk_log.txt
mv ${prefix}_LD_table_clean.txt ${prefix}_LD_table.txt

#tomahawk decay -i ${prefix}_LD.two -I 1:1-60000 > ${prefix}_decay.txt

####################################################
### 	    			 END 			 		 ###
####################################################
