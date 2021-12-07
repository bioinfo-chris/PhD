source ~/.zshrc
##################################################################
###   Pairwise SNP Linkage Disequilibrium Satistics From MSA   ###
##################################################################

# this script requires perl, snp-sites, bcftools and tomahawk to be installed and symlinked
# it will produce biallelic haploid and diploid .vcf files from a multiple sequence alignment (msa)
# as well as a diploid .bcf and tomahawk intermediate files (.twk, .two) and the final pairwise ld table (.txt)

# usage:
# bash msa2ld.sh <(path/to/)msa.fasta> <(output/path/)prefix> <ref_seq_header>

# get path/msa.fasta as first argument
input=$1

# get output prefix (and path, if specifying different output directory)
# if ignored, defaults to msa prefix and directory
prefix=$2
if [[ -z $prefix ]] ; then
	prefix=$(echo ${input} | rev | awk 'BEGIN{FS="."}{for (i=2; i<NF; i++) printf $i "."; print $NF}' | rev)
fi

# reference sequence header (same as in alignment)
# if ignored, default uses first sequence in msa as the reference
ref=$3

# set output names
hap_out=$(echo $prefix | sed 's/$/_haploid.vcf/')
dip_out=$(echo $prefix | sed 's/$/_diploid.vcf/')
bcf_out=$(echo $prefix | sed 's/$/_diploid.bcf/')

# make snp-sites vcf
if [[ -z $ref ]] ; then
	snp-sites -v -o ${hap_out} ${input}
else
	ref_aln_out=$(echo $input | sed 's/.fasta/_reorder.fasta/')
	ref_seq=$(fasta2tbl ${input} | grep "\b${ref}\b" | tbl2fasta | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d')
	msa_no_ref=$(fasta2tbl ${input} | grep -v "\b${ref}\b" | tbl2fasta | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d')
	echo ${ref_seq}'\n'${msa_no_ref} > ${ref_aln_out}
	snp-sites -v -o ${hap_out} ${ref_aln_out}
	#rm ${ref_aln_out}
fi

# damien's perl code to convert tri+ allelic sites and indels, and diploidize alleles from haploid vcf
# genotype column must contain only one character (default output of snp-sites)
# kept alt allele is the one most represented among isolates, in case of equality, alt chosen alphabeticaly: A C G T
perl -F'\t' -ane 'if(/^##/){print}elsif(/^#/){print}else{ $a++;
$F[6] = "."; $F[-1] =~ s/\r?\n//; 
$star = 0 ;
@alleles = split /,/,$F[4] ; 
foreach $j (0..$#alleles){ $alleles[$j] =~ s/\r?\n//; if( $alleles[$j] =~ /\*/ ){ $star = $j; $star++  }  }

$nb_alt_isolates = 0;

undef %alt_allele_count;
foreach $i (9..$#F){  if($F[$i] == $star ){ $F[$i] = 0 } ; if($F[$i] != 0){ $alt_allele_count{$F[$i]}++ } };

$alt_to_keep = 0;
foreach $key (sort { $alt_allele_count{$b} <=> $alt_allele_count{$a} or $a cmp $b } keys %alt_allele_count){
$alt_to_keep = $key;
print STDERR "Log : Alternate allele to keep is $alleles[($alt_to_keep-1)] at position $F[1]\n";
last;
} 

foreach $i (9..$#F){ 
if(length($F[$i]) > 1){ print STDERR "Fatal error 1, do not use output file : more than one character is used to define the genotype of strains at position $F[1] , script will not work.\n" }
 if($F[$i] != $alt_to_keep ){ $F[$i] = 0 } ; if($F[$i] != 0){ $nb_alt_isolates++ ; $F[$i] = 1 } ; $F[$i] = $F[$i] . "/" . $F[$i] };
# i now need to update the allele considered as unique ALT in the 5th column
$F[4] = $alleles[($alt_to_keep-1)] ;
$top = join("\t",@F[0..$#F]); print $top . "\n";
}' ${hap_out} > ${dip_out}

# make a cut down vcf as SNP pres/abs table
tail -n +4 ${dip_out} | sed 's/#//;s@0/0@0@g;s@1/1@1@g' > ${prefix}_SNP_pres-abs_table.txt


# convert diploid vcf to binary
bcftools view ${dip_out} -O b -o ${bcf_out}

# run tomahawk
tomahawk import -f -i ${bcf_out} -o ${prefix}_LD
tomahawk calc -u -i ${prefix}_LD.twk -o ${prefix}_LD -t 6 -r 0
tomahawk view -i ${prefix}_LD.two > ${prefix}_LD_table.txt
tail -n +14 ${prefix}_LD_table.txt > ${prefix}_LD_table_clean.txt
head -13 ${prefix}_LD_table.txt > ${prefix}_tomahawk_log.txt
mv ${prefix}_LD_table_clean.txt ${prefix}_LD_table.txt

#tomahawk decay -i ${prefix}_LD.two -I 1:1-60000 > ${prefix}_decay.txt

####################################################
### 	    		END 			 ###
####################################################
