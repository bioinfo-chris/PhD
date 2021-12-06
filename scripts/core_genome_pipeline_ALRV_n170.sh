##############################################################################
###  	  	  ALRV 170 Isolates Core Delineation/Extraction	           ###
##############################################################################
source ~/.zshrc

# start with prokka, put whole genome set in prokka directory
prokka=/Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/prokka
cd ${prokka} && mcd alrv_n170
cp $PWD/alrv_n170_alignable_seqs_seq_id.fasta ./

# get single seqs
mcd single_seqs
cat ../*.fasta | awk 'BEGIN{RS=">"}{filename=($1".fasta"); print ">" $0 > filename; close(filename)}'

# run prokka
up && mkdir prokka
for i in single_seqs/* ; do
	name=$(echo $i | sed 's@single_seqs/@@;s/.fasta//')
	mkdir prokka/${name}
	prokka --kingdom Viruses --prefix ${name} --outdir prokka/${name} ${i} --force --locustag ${name} >> prokka_log.txt 2>&1 
done

# change to roary, set up gff dir
roary=/Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/roary/alrv_roary
cd ${roary} && mkdir gff_n170
cp ${prokka}/alrv_n170/prokka/**/*.gff gff_n170/

# run roary
roary -e --mafft -p 10 -i 80 -v -f alrv_n170_80_serv/ gff_n170/*.gff

# scoary
mcd alrv_n170_80_serv/scoary
tree=/Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/raxml/n49_ORFs_blastn/alrv_n170/RAxML_bipartitions.alrv_n170_CGS_n49_ml-bs_REROOT.tree
scoary -g ../gene_presence_absence.csv -t alrv_n170_clades.csv -n ${tree} --collapse --no_pairwise

# make a core_genomes directory in roary project and set as home
mcd ../alrv_n170_80_serv/core_genomes
home_dir=/Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/roary/alrv_roary/alrv_n170_80_serv/core_genomes

# make key dirs
mkdir blast_db
mkdir core_tblastn
mkdir alignments

# put whole genomes in blast_db
cp $PWD/alrv_n170_alignable_seqs_seq_id.fasta blast_db/
makeblastdb -in blast_db/*.fasta -out blast_db/alrv_n170_db -dbtype nucl

# move into core_tblastn
cd core_blastn

# *** edit roary pan-genome reference sets to only contain core genes, put in core_tblastn *** #

# begin extracting core via tblastn ...
# make core gene directories
mcd single_orf_seqs
cat ../*.fasta | awk 'BEGIN{RS=">"}{filename=($1".fasta"); print ">" $0 > filename; close(filename)}'
for i in * ; do
	seq_name=$(echo $i | awk 'BEGIN{FS="."}{print $1}')
	mkdir $seq_name
	mv $i ${seq_name}/${seq_name}.fasta
done

#tblastn each orf seq against each WGS
blast_db=${home_dir}/blast_db/alrv_n170_db
WGS=${home_dir}/blast_db/alrv_n170_alignable_seqs_seq_id.fasta

for i in * ; do
	mkdir ${i}/individ_seqs
	blastn -db $blast_db \
		-query ${i}/${i}.fasta \
		-qcov_hsp_perc 80 \
        	-perc_identity 80 \
		-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
		-out ${i}/${i}_blast_out.txt
	cat ${i}/${i}_blast_out.txt | awk -v var="$i" 'BEGIN {RS="\n";FS=" "} {print $2"."var, $13}' | tbl2fasta \
		| sed -e 's/\(^>.*$\)/#\1#/' $1 | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' \
		> ${i}/${i}_isolate_seqs.fasta
	cat ${i}/${i}_isolate_seqs.fasta | awk -v var="$i" 'BEGIN{RS=">"}{filename=(var"/individ_seqs/"$1".fasta"); print ">" $0 > filename; close(filename)}'
done

# count seqs
for i in * ; do
	echo "$i\t" $(ls ${i}/individ_seqs | wc -l) >> ../counts.txt
done
up

# make all seqs, cat and unify directories
mkdir all_complete_genes && mkdir seqs_catted_core_genes && mkdir seqs_unified_core_genes

# *** find genes with all isolates from counts.txt file, make list orfs_all_isos.txt ***
while read orf ; do
	 cp single_orf_seqs/${orf}/individ_seqs/* all_complete_genes/
done < orfs_all_isos.txt

# make ids file, cat all genes by isolate
grep ">" $WGS | sed 's/>//' > ids.txt
while read id ; do
		cat all_complete_genes/*${id}.* > seqs_catted_core_genes/${id}_catted_core.fasta
		union -filter seqs_catted_core_genes/${id}_catted_core.fasta > seqs_unified_core_genes/${id}_unified_core.fasta
		sed -i '' -e 's/.BALF5//' seqs_unified_core_genes/${id}_unified_core.fasta
done < ids.txt 
up

# move final CGS set and prepare for alignment
n_orfs=$(cat core_blastn/orfs_all_isos.txt | wc -l | awk '{print $1}')
mkdir alignments/n${n_orfs}_ORFs
cat core_blastn/seqs_unified_core_genes/* > alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs.fasta
stripn alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs.fasta > alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs.fasta_
mv alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs.fasta_ alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs.fasta

# align and trim
mafft --thread 10 --reorder alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs.fasta > alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs_aln.fasta
trimal -in alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs_aln.fasta -out alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20.fasta -gt 0.8

# build raxml tree
raxml=/Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/raxml
mkdir -p ${raxml}/n${n_orfs}_ORFs_blastn/alrv_n170
cd ${raxml}/n${n_orfs}_ORFs_blastn/alrv_n170

# first find & filter homoplasies
mkdir mpboot
mkdir homoplasy_finder

# set up and run mpboot
cp ${home_dir}/alignments/n${n_orfs}_ORFs/alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20.fasta ./
cp *.fasta mpboot/
mpboot -s mpboot/*.fasta -bb 1000
rm mpboot/*.fasta

# run homoplasy finder
homofinder=/Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/scripts/homoplasy_screen_script.R
Rscript --vanilla ${homofinder} -m mpboot/*.treefile -a *.fasta -o homoplasy_finder/
cp homoplasy_finder/*.fasta ./alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20_nohomo.fasta
snp-sites -o alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20_SNPs.fasta alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20.fasta
snp-sites -o alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20_nohomo_SNPs.fasta alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20_nohomo.fasta

# edit homoplasy outputs
mv homoplasy_finder/*consistencyIndexReport* homoplasy_finder/consistencyIndexReport.txt
echo 'alrv homoplasies: '$(( $(wc -l homoplasy_finder/*.txt | awk '{print $1}') - 1)) > homoplasy_finder/homoplasy_count.txt
echo 'alrv SNPs before: '$(seqlen alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20_SNPs.fasta | head -1 | awk '{print $2}') >> homoplasy_finder/homoplasy_count.txt
echo 'alrv SNPs after: '$(seqlen alrv_n170_CGS_n${n_orfs}_ORFs_aln_trim20_nohomo_SNPs.fasta | head -1 | awk '{print $2}') >> homoplasy_finder/homoplasy_count.txt

# set up and run raxml
mkdir bs_tree && mcd ml_tree 
raxmlHPC-SSE3 -m GTRCAT -s ../*nohomo.fasta -n alrv_n170_CGS_n${n_orfs}_nohomo_ml_tree -# 100 -p 12345 -T 20
cd ../bs_tree
raxmlHPC-SSE3 -m GTRCAT -s ../*nohomo.fasta -n alrv_n170_CGS_n${n_orfs}_nohomo_bs_tree -# 1000 -p 12345 -b 12345 -T 20 
up
raxmlHPC-SSE3 -f b -m GTRCAT -z bs_tree/RAxML_bootstrap.* -t ml_tree/RAxML_bestTree.* -n alrv_n170_CGS_n${n_orfs}_ml-bs.tree


##############################################################################
###  	  	   		     END	       			   ###
##############################################################################
