#!/bin/bash
# set pointer to program.
program="silp.py"
work_dir=".."
ref_dir="${work_dir}/ref"
asm_dir="${work_dir}/asm"
aln_dir="${work_dir}/aln"
scf_dir="${work_dir}/scf"

# align
python $program align \
	-w $scf_dir \
	-a $aln_dir \
	-p 5 \
	-c ${asm_dir}/asm.fasta \
	-s ${asm_dir}/asm.length \
	-q1 ${ref_dir}/read1.fastq \
	-q2 ${ref_dir}/read2.fastq \
	-k 2

# preprocess
python $program nodes -w ${scf_dir} -c ${asm_dir}/asm.fasta
python $program edges -w ${scf_dir} -i 3500 -s 350 -rf -s1 ${aln_dir}/read1.sam -s2 ${aln_dir}/read2.sam
python $program bundles -w ${scf_dir} -b 2 -p 90 -bup 1 -r ${aln_dir}/ant -i 3500 -s 350
python $program decompose -w ${scf_dir} -m 2500

# start time
start=$(date +%s)

# scaffold
python $program orient -w $scf_dir -z 0
python $program order -w $scf_dir
python $program gap -w $scf_dir
python $program write -w $scf_dir -a ${scf_dir}/scf.agp
python $program fasta -w $scf_dir -a ${scf_dir}/scf.agp -c ${asm_dir}/asm.fasta -f ${scf_dir}/scf.fasta

# stop time
stop=$(date +%s)
echo RUNTIME: $(expr $stop - $start) >> ${scf_dir}/scf.agp

