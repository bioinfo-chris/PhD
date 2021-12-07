###############################################################
###		 Linkage Disequilibrium Pipeline            ###
###############################################################
# call script as zsh *not* bash
source ~/.zshrc
zmodload zsh/zutil

# this script reconstructs a set of core genomes (one multifasta per core genome) so that ORFs are in same order and direction of specified reference whole genome(s)
# it then runs a series of scripts to perform and plot LD calculations, which must be in same directory as this script
# these are: msa2ld.sh ; vcf2ld.sh ; var_calc.R ; ld_perms_and_plots.R

# usage:
# zsh LD_pipeline.sh 		-h --home-dir       <str> (OPTIONAL)	specify a home directory, must contain additional scripts
#									if not provided, directory where script is used
#				-p --project        <str> (REQUIRED)	specific system project directory name which will be created in home directory 
#				-i --id-list        <str> (REQUIRED)	path to list of sequence/isolate ids (without .fasta suffix) matching input wgs and core, used to fetch them 
#				-w --wgs-input      <str> (REQUIRED)	path to directory of all single whole genome sequence assemblies, fetch using id list
# 									e.g. /Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/genomes/complete_genomes/all_whole_genomes/final_seqs/single_seqs
#				-c --core-input     <str> (REQUIRED)	path to directory of all single catted (not unified) core genome sequences 
#									must match wgs names with added suffix <.core_genes>
#									e.g. /Users/ChrisOwen/Dropbox/Work/PhD/Chapters/3-Ranavirus/roary/atv_roary/atv_80_ehnv/core/isolates/isolate_core_genes
# 				-r --reference-list <str> (OPTIONAL)	path to list of inpust sequences to be used as references
#									if not supplied, useses all input sequence as references
#				-f --freq-list 	    <str> (OPTIONAL)	path to list of minimum threshold frequencies of varient sites to filter out
#									if not supplied, will filter sits under 0.03, 0.04, 0.08, 0.10, 0.16

# get arguments
zparseopts -D -E - h:=home_dir -home_dir:=home_dir \
	p:=proj_dir -project:=proj_dir \
	w:=input_wgs -wgs-input:=input_wgs \
	c:=input_core -core-input:=input_core \
	i:=ids -id-list:=ids \
	r:=ref_list -reference-list:=ref_list \
	f:=var_freqs -freq-list:=var_freqs || exit 1

# make arguments as should be
home_dir=${home_dir[-1]}
proj_dir=${proj_dir[-1]}
input_wgs=${input_wgs[-1]}
input_core=${input_core[-1]}
ref_list=${ref_list[-1]}
ids=${ids[-1]}
var_freqs=${var_freqs[-1]}

# BEGIN
echo "=============================================="
echo "      ----===::: LD Pipeline :::===----       "
echo "=============================================="
echo ""
date
echo ""

# change to home directory, if not provided, home is where script is
# (ideally create LD analysis home directory and put/run all scripts from there)
if [[ -z ${home_dir} ]] ; then
	home_dir=$PWD
else
	mkdir ${home_dir}
	cd ${home_dir}
fi

# make and cd into project directory, get ids
proj_check=$(ls | grep "\\\b${proj_dir}\\\b") 
if [[ -z ${proj_check} ]] ; then
	mkdir ${proj_dir} && cd ${proj_dir}
	cat ${ids} > ids.txt
else
	cd ${proj_dir}
	cat ${ids} > ids.txt
fi

# check for var freqs option
if [[ -z $var_freqs ]] ; then
	echo 0.03'\n'0.04'\n'0.08'\n'0.10'\n'0.16 > var_freqs.txt
else
	cat ${var_freqs} > var_freqs.txt
fi

# check if the key dirs (and therefore seqs) are present
# if not, make key sub dirs
# copy wgs and their core to relevent directory based on list of ids
seq_check=$(ls | grep "core")
if [[ -z ${seq_check} ]] ; then
	echo "making key project directories and copying genomes ..." 
	mkdir core && mkdir ref_alns && mkdir WGS
	mkdir -p core/isolate_core_genes && mkdir -p core/all_genes_all_seqs
	while read iso ; do
		cp ${input_wgs}/${iso}.fasta WGS/
		cp ${input_core}/${iso}.core_genes.fasta core/isolate_core_genes/
	done < ids.txt
fi

# split into individ genes
seq_check=$(ls core/all_genes_all_seqs/)
if [[ -z ${seq_check} ]] ; then 
	cd core/all_genes_all_seqs
	for i in ../isolate_core_genes/* ; do
		cat $i | awk 'BEGIN{RS=">"}{filename=($1".fasta"); print ">" $0 > filename; close(filename)}'
	done
	cd ${home_dir}/${proj_dir}
	echo "done"
	echo ""
fi

# make list of ref dirs to make
if [[ -z ${ref_list} ]] ; then
	cat ${home_dir}/ids.txt | sed 's/.fasta//;s/^/ref_/' > ref_list.txt
else
	cat ${ref_list} > ref_list.txt
fi
it_tot=$(cat ref_list.txt | wc -l | awk '{print $1}')

# get going with refs ...
cd ref_alns

# make ref project directories
it=0 ; while read ref ; do
	it=$((${it}+1))
	ref_head=$(echo ${ref} | sed 's/ref_//')
	echo "=============================================="
	echo "processing isolates according reference ${ref_head} (${it} of ${it_tot}) ..."
	echo ""
	echo "setting up reference directories ..."
	init_check=$(ls ${ref}/ | grep "tomahawk") 2> /dev/null
	if [[ -z ${init_check} ]] ; then
		mkdir -p ${ref}/alignments
		mkdir ${ref}/blast_check
		mkdir ${ref}/WGS_blast
		mkdir -p ${ref}/tomahawk/all_vars
		mkdir -p ${ref}/core_reorder/rev_comp
		mkdir ${ref}/core_reorder/isolates
	fi
	echo "done"

	# perform wgs blast against all its core
	init_check=$(ls ${ref}/WGS_blast |  grep "_blast_out.txt")
	if [[ -z ${init_check} ]] ; then
		makeblastdb -in ${home_dir}/${proj_dir}/core/isolate_core_genes/${ref_head}.core_genes.fasta -out ${ref}/WGS_blast/${ref_head}_db -dbtype nucl > /dev/null 2>&1 
		blastn -db ${ref}/WGS_blast/${ref_head}_db \
				-query ${home_dir}/${proj_dir}/WGS/${ref_head}.fasta \
			    	-perc_identity 90 \
				-outfmt 6 \
				-out ${ref}/WGS_blast/${ref_head}_blast_out.txt
		cat ${ref}/WGS_blast/${ref_head}_blast_out.txt | sort -rnk4 | sort -uk2,2 | column -t > ${ref}/WGS_blast/${ref_head}_blast_out.txt_
		mv ${ref}/WGS_blast/${ref_head}_blast_out.txt_ ${ref}/WGS_blast/${ref_head}_blast_out.txt
	fi

	# based on blasts of ref against ref core, make list of reverse compliment genes realtive to ref wgs
	# use list to rev comp genes that need it
	echo "sorting out genes in relation to reference ..."
	init_check=$(ls ${ref}/core_reorder |  grep "all_isolates_core_${ref}_reordered.fasta")
	if [[ -z ${init_check} ]] ; then
		while read blast ; do
			gene=$(echo $blast | awk '{print $2}' | rev | awk 'BEGIN{FS="."}{print $1}' | rev)
			ori=$(echo $blast | awk '{print $10}')
			if [ $ori = "1" ] ; then
				echo $gene | sed 's/^/\\\b/;s/$/\\\b/' >> ${ref}/core_reorder/${ref}_rev_comp_list.txt
			fi
		done < ${ref}/WGS_blast/*_blast_out.txt
		ls ${home_dir}/${proj_dir}/core/all_genes_all_seqs | grep -f ${ref}/core_reorder/${ref}_rev_comp_list.txt | sed 's/.fasta//' > ${ref}/core_reorder/${ref}_seqs_to_rev.txt
		it2=0
		it_tot2=$(cat ${ref}/core_reorder/${ref}_seqs_to_rev.txt | wc -l | awk '{print $1}')
		while read rev ; do
			it2=$((${it2}+1))
			perc2=$((100*${it2}/${it_tot2}))
			if [ ${it2} -ne ${it_tot2} ] ; then
				echo -ne " reverse complimenting sequences ... (${perc2}%) \r"
			else
				echo -ne " reverse complimenting sequences ... (99%) \r"
			fi
			revseq -sequence ${home_dir}/${proj_dir}/core/all_genes_all_seqs/${rev}.fasta -outseq ${ref}/core_reorder/rev_comp/${rev}.revcomp.fasta -notag 2> /dev/null
			stripn ${ref}/core_reorder/rev_comp/${rev}.revcomp.fasta > ${ref}/core_reorder/rev_comp/${rev}.revcomp.fasta_
			mv ${ref}/core_reorder/rev_comp/${rev}.revcomp.fasta_ ${ref}/core_reorder/rev_comp/${rev}.revcomp.fasta
			if [ ${it2} -eq ${it_tot2} ] ; then
				echo " reverse complimenting sequences ... (100%)"
			fi
		done < ${ref}/core_reorder/${ref}_seqs_to_rev.txt
		cp -r ${home_dir}/${proj_dir}/core/all_genes_all_seqs/ ${ref}/core_reorder/all_genes_for/
		while read seq ; do
			rm ${ref}/core_reorder/all_genes_for/${seq}.fasta
		done < ${ref}/core_reorder/${ref}_seqs_to_rev.txt
		cp ${ref}/core_reorder/rev_comp/* ${ref}/core_reorder/all_genes_for/
	fi
	echo "done"

	# copy forward genes to isolate specific directories
	init_check=$(ls ${ref}/core_reorder |  grep "all_isolates_core_${ref}_reordered.fasta")
	if [[ -z ${init_check} ]] ; then
		while read iso ; do
			mkdir -p ${ref}/core_reorder/isolates/${iso}/single_gene_seqs
			cp ${ref}/core_reorder/all_genes_for/*${iso}.* ${ref}/core_reorder/isolates/${iso}/single_gene_seqs/
		done < ../ids.txt
	fi

	# get gene order list from blast out
	init_check=$(ls ${ref}/alignments |  grep "all_isolates_core_${ref}_reordered_WGS_positions_aln.fasta")
	if [[ -z ${init_check} ]] ; then
		blast_sort=$(cat ${ref}/WGS_blast/*_blast_out.txt | sort -nk 7)
		echo $blast_sort | awk '{print $2}' | rev | awk 'BEGIN{FS="."}{print $1}' | rev > ${ref}/core_reorder/${ref}_gene_order.txt
	fi

	# concatenate core genes in order that they occour in ref, then unify
	echo "building isolate core genomes according to reference ..."
	init_check=$(ls ${ref}/core_reorder |  grep "all_isolates_core_${ref}_reordered.fasta")
	if [[ -z ${init_check} ]] ; then
		it2=0
		it_tot2=$(ls ${ref}/core_reorder/isolates/ | wc -l | awk '{print $1}')
		for i in ${ref}/core_reorder/isolates/* ; do
			it2=$((${it2}+1))
			perc2=$((100*${it2}/${it_tot2}))
			if [ ${it2} -ne ${it_tot2} ] ; then
				echo -ne " concatenating and unifying core genes ... (${perc2}%) \r"
			else
				echo -ne " concatenating and unifying core genes ... (99%) \r"
			fi
			iso=$(echo $i | sed "s@${ref}/core_reorder/isolates/@@")
			while read order ; do
				cat ${i}/single_gene_seqs/*${order}.* >> ${i}/${iso}.core_ref_order-direc.fasta
			done < ${ref}/core_reorder/${ref}_gene_order.txt
			union -filter ${i}/${iso}.core_ref_order-direc.fasta > ${i}/${iso}.core_ref_order-direc_unified.fasta
			sed -i '' -e "1 s/^.*$/>${iso}/" $i/${iso}.core_ref_order-direc_unified.fasta
			stripn ${i}/${iso}.core_ref_order-direc_unified.fasta > ${i}/${iso}.core_ref_order-direc_unified.fasta_
			mv ${i}/${iso}.core_ref_order-direc_unified.fasta_ ${i}/${iso}.core_ref_order-direc_unified.fasta
			if [ ${it2} -eq ${it_tot2} ] ; then
				echo " concatenating and unifying core genes ... (100%)"
			fi
		done
	fi
	echo "done"

	# cat into multifasta of all core genomes
	init_check=$(ls ${ref}/core_reorder |  grep "all_isolates_core_${ref}_reordered.fasta")
	if [[ -z ${init_check} ]] ; then
		for iso in ${ref}/core_reorder/isolates/* ; do
			cat ${iso}/*unified.fasta >> ${ref}/core_reorder/all_isolates_core_${ref}_reordered.fasta
		done
	fi

	# align core, trim, align to ref wgs, remove wgs
	echo "aligning isolate core genomes in whole genome positions  ..."
	init_check=$(ls ${ref}/alignments |  grep "all_isolates_core_${ref}_reordered_WGS_positions_aln.fasta")
	if [[ -z ${init_check} ]] ; then
		mafft ${ref}/core_reorder/all_isolates_core_${ref}_reordered.fasta > ${ref}/alignments/all_isolates_core_${ref}_reordered_aln.fasta 2> ${ref}/alignments/${ref}_aln_log.txt
		trimal -in ${ref}/alignments/all_isolates_core_${ref}_reordered_aln.fasta -out ${ref}/alignments/all_isolates_core_${ref}_reordered_aln_trim20.fasta -gt 0.8 >> ${ref}/alignments/${ref}_aln_log.txt 2>&1
		wgs=${home_dir}/${proj_dir}/WGS/${ref_head}.fasta
		muscle -profile -in1 ${wgs} -in2 ${ref}/alignments/all_isolates_core_${ref}_reordered_aln_trim20.fasta -out ${ref}/alignments/all_isolates_core_${ref}_reordered_with-WGS_aln.fasta >> ${ref}/alignments/${ref}_aln_log.txt 2>&1
		stripn ${ref}/alignments/all_isolates_core_${ref}_reordered_with-WGS_aln.fasta > ${ref}/alignments/all_isolates_core_${ref}_reordered_with-WGS_aln.fasta_
		mv ${ref}/alignments/all_isolates_core_${ref}_reordered_with-WGS_aln.fasta_ ${ref}/alignments/all_isolates_core_${ref}_reordered_with-WGS_aln.fasta
		tail -n +3 ${ref}/alignments/all_isolates_core_${ref}_reordered_with-WGS_aln.fasta > ${ref}/alignments/all_isolates_core_${ref}_reordered_WGS_positions_aln.fasta
	fi
	echo "done"

	# run a blast of ref wgs against the ref core genes, check order and position
	init_check=$(ls ${ref}/blast_check |  grep "_blast_check_out.txt")
	if [[ -z ${init_check} ]] ; then
		makeblastdb -in ${ref}/core_reorder/isolates/${ref_head}/${ref_head}.core_ref_order-direc.fasta -out ${ref}/blast_check/${ref_head}_db -dbtype nucl > /dev/null 2>&1 
		blastn -db ${ref}/blast_check/${ref_head}_db \
					-query ${home_dir}/${proj_dir}/WGS/${ref_head}.fasta \
			        	-perc_identity 90 \
					-outfmt 6 \
					-out ${ref}/blast_check/${ref_head}_blast_check_out.txt
		cat ${ref}/blast_check/${ref_head}_blast_check_out.txt | sort -rnk4 | sort -uk2,2 | column -t > ${ref}/blast_check/${ref_head}_blast_check_out.txt_
		mv ${ref}/blast_check/${ref_head}_blast_check_out.txt_ ${ref}/blast_check/${ref_head}_blast_check_out.txt
	fi

	# use msa2ld.sh script to run LD analysis, R script to calculate varient freqs, filter vcf, re-run tomahawk on subsets, and plot results
	echo "calculating linkage disequilibrium statistics ..."
	init_check=$(ls ${ref}/tomahawk/ | grep -F "all_vars_n")
	if [[ -z ${init_check} ]] ; then
		ref_head=$(echo ${ref} | sed 's/ref_//')
		zsh ${home_dir}/msa2ld.sh ${ref}/alignments/all_isolates_core_${ref}_reordered_WGS_positions_aln.fasta ${ref}/tomahawk/all_vars/all_isolates_core_${ref}_reordered_WGS_positions ${ref_head} > ${ref}/tomahawk/all_vars/${ref}_msa2ld_log.txt 2>&1
		count=$(tail -n +2 ${ref}/tomahawk/all_vars/*_SNP_pres-abs_table.txt | wc -l | awk '{print $1}')
		mv ${ref}/tomahawk/all_vars ${ref}/tomahawk/all_vars_n${count}
		Rscript --vanilla ${home_dir}/var_calc.R ${proj_dir}/ref_alns/${ref} ${home_dir}/${proj_dir}/var_freqs.txt
		sed -i '' -e 's/^/\\\b/;s/$/\\\b/' ${ref}/tomahawk/vars_*/*.txt
		for freq in ${ref}/tomahawk/vars_gt* ; do
			freqn=$(echo $freq | sed "s@${ref}/tomahawk/vars_@@")
			head -4 ${ref}/tomahawk/all_vars_n*/*_diploid.vcf >> ${freq}/all_isolates_core_${ref}_reordered_WGS_positions_${freqn}_diploid.vcf
			tail -n +5 ${ref}/tomahawk/all_vars_n*/*_diploid.vcf | grep -f ${freq}/*.txt >> ${freq}/all_isolates_core_${ref}_reordered_WGS_positions_${freqn}_diploid.vcf
			zsh ${home_dir}/vcf2ld.sh ${freq}/all_isolates_core_${ref}_reordered_WGS_positions_${freqn}_diploid.vcf > ${freq}/${ref}_${freqn}_vcf2ld_log.txt 2>&1
		done
		Rscript --vanilla ${home_dir}/ld_perms_and_plots.R ${proj_dir}/ref_alns/${ref}
	fi
	echo "done"
	echo ""
	echo "finished analysis according to reference ${ref} (${it} of ${it_tot})"
	echo "=============================================="
done < ../ref_list.txt
echo ""
echo "all references done"
echo ""

# harvest the r2 regression R-squared and coefficient results
echo "ref_seq\tvars\tperm_pvalue\tcoef\tr2" > ${home_dir}/${proj_dir}/ref_r2_regression_all_coefs-rsq.txt
for ref in * ; do
	for freq_dir in ${ref}/tomahawk/vars_gt* ; do
		freqn=$(echo $freq_dir | sed "s@${ref}/tomahawk/vars_gt@@")
		freq=$(echo $freqn | awk 'BEGIN{FS="_n"}{print $1}')
		pvalue=$(cat ${ref}/tomahawk/vars_gt${freqn}*/*_permutation_pvalue.txt | awk 'BEGIN{FS=": "}{print $2}')
		r2=$(cat ${ref}/tomahawk/vars_gt${freqn}*/*r2-dist_lm.txt | tail -3 | head -1 | awk 'BEGIN{FS=":"}{print $3}' | awk '{print $1}')
		coef=$(cat ${ref}/tomahawk/vars_gt${freqn}*/*r2-dist_lm.txt | tail -8 | head -1 | awk 'BEGIN{FS="dist.pos"}{print $2}' | awk '{print $1}')
		echo ${ref}'\t'vars_gt${freq}'\t'${pvalue}'\t'${coef}'\t'${r2} >> ${home_dir}/${proj_dir}/ref_r2_regression_all_coefs-rsq.txt
	done
done
cat ${home_dir}/${proj_dir}/ref_r2_regression_all_coefs-rsq.txt | column -t | sed $'s/[^[:print:]\t]//g' > ${home_dir}/${proj_dir}/ref_r2_regression_all_coefs-rsq.txt_
mv ${home_dir}/${proj_dir}/ref_r2_regression_all_coefs-rsq.txt_ ${home_dir}/${proj_dir}/ref_r2_regression_all_coefs-rsq.txt

# clean up
echo "cleaning up ..."
it=0 ; while read ref ; do
	init_check=$(ls ${ref}/core_reorder/ | grep "all_genes_for")
	if [[ ${init_check} ]] ; then
		it=$((${it}+1))
		perc=$((100*${it}/${it_tot}))
		ref_head=$(echo ${ref} | sed 's/ref_//')
		if [ ${it} -ne ${it_tot} ] ; then
			echo -ne "removing trash from reference ${ref_head} (${it} of ${it_tot}) ... (${perc}%) \r"
		else
			echo -ne "removing trash from reference ${ref_head} (${it} of ${it_tot}) ... (99%) \r"
		fi
		rm -r ${ref}/core_reorder/all_genes_for
		rm -r ${ref}/core_reorder/rev_comp
		rm ${ref}/tomahawk/**/*.two ${ref}/tomahawk/**/*.twk
		if [ ${it} -eq ${it_tot} ] ; then
			echo "removing trash from reference ${ref_head} (${it} of ${it_tot}) ... (100%)"
		fi
	fi
done < ../ref_list.txt
echo ""
echo "done & done"
echo ""
date

###############################################################
###			     END			    ###
###############################################################
