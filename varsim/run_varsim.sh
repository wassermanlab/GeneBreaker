#!/bin/bash

## Account Information
#SBATCH --account=rrg-wyeth
#SBATCH --job-name=varsim_testrun

## Mail Options
#SBATCH --mail-user=avshalom.tamar0@gmail.com
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.output
#SBATCH --error=%x-%j.error


module load java/1.8.0_121
source activate variant_sim
cd /project/projects/def-wyeth/avshalom/variant_simulation_framework/varsim_run

conda=`conda info |grep 'envs directories' |awk '{print $4}'`
bin="${conda}/variant_sim/bin/"

./varsim.py --vc_in_vcf All.vcf.gz --sv_insert_seq insert_seq.txt \
--sv_dgv GRCh37_hg19_supportingvariants_2013-07-23.txt \
--reference hs37d5.fa --id simu --read_length 150 --vc_num_snp 0 --vc_num_ins 0 \
--vc_num_del 0 --vc_num_mnp 0 --vc_num_complex 0 --sv_num_ins 0 \
--sv_num_del 0 --sv_num_dup 0 --sv_num_inv 0 --sv_percent_novel 0 \
--vc_percent_novel 0 --mean_fragment_size 450 --sd_fragment_size 50 \
--vc_min_length_lim 0 --vc_max_length_lim 0 --sv_min_length_lim 0 \
--sv_max_length_lim 0 --nlanes 1 --total_coverage 30 \
--simulator_executable ${bin}art_illumina --out_dir out --log_dir log --work_dir work \
--simulator art \
--vcfs /project/projects/def-wyeth/FastTyper/simulations/vcfs/0simulation.vcf \
--art_options "-ef"
