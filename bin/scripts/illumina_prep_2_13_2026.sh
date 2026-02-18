#!/bin/bash
#SBATCH --time=48:00:00
proj_env=$1
NAME=$2
Illumina_meta_csv=$3
source $proj_env

out_dir="${my_results}/${NAME}"
prep_dir="${out_dir}/Illumina"
mkdir -p ${out_dir} ${prep_dir}
echo 'Making Directories for New Genotypes'
echo 'Copying Genotypes from Source Directory'


look for paths in Illumina_meta_csv see xample below
full_path
/blue/mateescu/raluca/UF Genotypes Feb2024/Univ_of_Florida_Marsella_BOVF250V1_20211108_FinalReport.txt
/blue/mateescu/raluca/UF Genotypes Feb2024/Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt
/blue/mateescu/raluca/UF Genotypes Feb2024/Univ_of_Florida_Mateescu_BOVF250V1_20160401_FinalReport.txt
/blue/mateescu/raluca/UF Genotypes Feb2024/Univ_of_Florida_Mateescu_BOVF250V1_20160407_FinalReport.txt
/blue/mateescu/raluca/UF Genotypes Feb2024/Univ_of_Florida_Mateescu_BOVF250V1_20160425_FinalReport.txt

cp them into the current dir


# Before we proceed we will check the raw Illumina files and take notes on:
# Now I recommend manually checking each header of each file and making sure they follow the same format of the example in the Protocol
head -20 ${prep_dir}/*.txt
## The number of samples per file and total number of samples
grep "Num Samples" Univ_of_Florida_*.txt
# This Grep function will find the total number of Samples per file
# The following for loop goes through all the Illumina Gentoype files and prints the number of columns in the files
# The first number is the smallest column number which should be one
# The 3rd number is the largest column number which should be 26 belonging to the header row of the Illumina genotype files
# The 2nd number, is the number of interest and should be 11. This corresponds to the columns in the actual genotypes portion of the file. If this number is not 11, this genotype file will have to be manipulated/reformatted.
cd /blue/mateescu/raluca/'UF Genotypes Feb2024'/
for file in `ls Univ_of_Florida_*BOVF250*.txt`
  do 
  echo ${file}
  head -20 ${file} |awk '{print NF}' | sort -nu | head -n 1
  head -20 ${file} |awk '{print NF}' | sort -nu | tail -n 2
  done
  
# At today's date (Jan 2022) the only file with this error is genotype file Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.tx. 
# This file contained 12 columns instead of 11. Due the Genotype ID column contatining an extra space. To correct this see the script below, which creates a temporary file in your new Illumina folder.
# This file combines columns 2 and 3 using a "_" and then replaces the file in the NEW Illumina folder.
# I do not recommend replacing the raw files from the source directory.
cd ${out_dir}/${Geno_name_date}/Illumina
#awk -F " " 'BEGIN { OFS="\t" } {print $1,$2_$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt > temp.txt

#mv temp.txt  Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt
#head -20 Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt
# Note this error is unique to this file and other files in the future may have different formatting issues. It is on you to figure out what the issue is and to check and fix that issue.
# Using this same line of code will not fix issues in future files. 
#Once all of this is done you can proceed to Illumina_to_plink.sh which contains code to convert these Illumina files into PLINK Files.s
#cd ${out_dir}/${Geno_name_date}/Illumina
#sed -i 's/Dc 14/Dc_14/g;s/Dc 116/Dc_116/g;s/BN 186/BN_186/g;s/Dc 501/Dc_501/g;s/Si 163/Si_163/g;s/Bn 34/Bn_34/g;s/Dc 125/Dc_125/g;s/Bn 175/Bn_175/g;s/Bn 33/Bn_33/g;s/Si 162/Si_162/g;s/Si 161/Si_161/g;s/Bn 164/Bn_164/g;s/Dc 529/Dc_529/g;s/Si 15/Si_15/g;s/Bn 157/Bn_157/g;s/Si 19/Si_19/g;s/Bn 173/Bn_173/g;s/Dc 528/Dc_528/g;s/Si 20/Si_20/g;s/Dc 128/Dc_128/g;s/Dc 543/Dc_543/g;s/Si 134/Si_134/g;s/Dc 16/Dc_16/g;s/Dc 15/Dc_15/g;s/Dc 518/Dc_518/g;s/DC 126/DC_126/g;s/SI 160/SI_160/g;s/DC 106/DC_106/g;s/SI 137/SI_137/g;s/DC 542/DC_542/g;s/Dc 12/Dc_12/g;s/BN 158/BN_158/g;s/Bn 179/Bn_179/g;s/DC 17/DC_17/g;s/SI 140/SI_140/g;s/SI 150/SI_150/g;s/DC 113/DC_113/g;s/DC 107/DC_107/g;s/SI 151/SI_151/g;s/DC 527/DC_527/g;s/DC 19/DC_19/g;s/DC 118/DC_118/g;s/SI 144/SI_144/g;s/Dc 120/Dc_120/g;s/Dc 509/Dc_509/g;s/Dc 524/Dc_524/g;s/Dc 122/Dc_122/g;s/Dc 121/Dc_121/g;s/Dc 20/Dc_20/g;s/Dc 112/Dc_112/g;s/Bn 23/Bn_23/g;s/Dc 549/Dc_549/g;s/Bn 32/Bn_32/g;s/DC 532/DC_532/g;s/Bn 35/Bn_35/g;s/Si 21/Si_21/g;s/Dc 551/Dc_551/g;s/Dc 103/Dc_103/g;s/Dc 550/Dc_550/g;s/Dc 105/Dc_105/g;s/Bn 27/Bn_27/g;s/Bn 24/Bn_24/g;s/Bn 174/Bn_174/g;s/Si 16/Si_16/g;s/Dc 500/Dc_500/g;s/Bn 153/Bn_153/g;s/Dc 130/Dc_130/g;s/Bn 629/Bn_629/g;s/Dc 109/Dc_109/g;s/Dc 29/Dc_29/g;s/Bn 185/Bn_185/g;s/Bn 152/Bn_152/g;s/Si 157/Si_157/g;s/Si 136/Si_136/g;s/Si 158/Si_158/g;s/Dc 18/Dc_18/g;s/Si 152/Si_152/g;s/Bn 159/Bn_159/g;s/Si 18/Si_18/g;s/Si 139/Si_139/g;s/Bn 168/Bn_168/g;s/Si 146/Si_146/g;s/Bn 31/Bn_31/g;s/Dc 25/Dc_25/g;s/Si 154/Si_154/g;s/Dc 510/Dc_510/g;s/Si 166/Si_166/g;s/Si 133/Si_133/g;s/Si 135/Si_135/g;s/Dc 124/Dc_124/g;s/Si 141/Si_141/g;s/Si 155/Si_155/g;s/SI 138/SI_138/g;s/BN 155/BN_155/g;s/SI 159/SI_159/g;s/SI 13/SI_13/g;s/SI 17/SI_17/g;s/BN 169/BN_169/g;s/BN 151/BN_151/g;s/BN 154/BN_154/g;s/SI 145/SI_145/g;s/DC 546/DC_546/g;s/SI 149/SI_149/g;s/SI 117/SI_117/g;s/DC 114/DC_114/g;s/DC 27/DC_27/g;s/BN 163/BN_163/g;s/SI 143/SI_143/g;s/DC 540/DC_540/g;s/BN 25/BN_25/g;s/SI 14/SI_14/g;s/54G C/54G_C/g;s/79IA C/79IA_C/g;s/47BR 2343/47BR_2343/g;s/301-9470/301-9470/g;s/DC119/DC119/g;s/8  8 Jan/8__8_Jan/g;s/BR 85/BR_85/g;s/UF 3/UF_3/g;s/BR 121/BR_121/g;s/UF 3/UF_3/g;s/79IA C/79IA_C/g' Univ_of_Florida_Mateescu_BOVF250V1_20220505_FinalReport.txt
#cd ${out_dir}/${Geno_name_date}/Illumina
#sed -i 's/Dc 14/Dc_14/g;s/Dc 116/Dc_116/g;s/BN 186/BN_186/g;s/Dc 501/Dc_501/g;s/Si 163/Si_163/g;s/Bn 34/Bn_34/g;s/Dc 125/Dc_125/g;s/Bn 175/Bn_175/g;s/Bn 33/Bn_33/g;s/Si 162/Si_162/g;s/Si 161/Si_161/g;s/Bn 164/Bn_164/g;s/Dc 529/Dc_529/g;s/Si 15/Si_15/g;s/Bn 157/Bn_157/g;s/Si 19/Si_19/g;s/Bn 173/Bn_173/g;s/Dc 528/Dc_528/g;s/Si 20/Si_20/g;s/Dc 128/Dc_128/g;s/Dc 543/Dc_543/g;s/Si 134/Si_134/g;s/Dc 16/Dc_16/g;s/Dc 15/Dc_15/g;s/Dc 518/Dc_518/g;s/DC 126/DC_126/g;s/SI 160/SI_160/g;s/DC 106/DC_106/g;s/SI 137/SI_137/g;s/DC 542/DC_542/g;s/Dc 12/Dc_12/g;s/BN 158/BN_158/g;s/Bn 179/Bn_179/g;s/DC 17/DC_17/g;s/SI 140/SI_140/g;s/SI 150/SI_150/g;s/DC 113/DC_113/g;s/DC 107/DC_107/g;s/SI 151/SI_151/g;s/DC 527/DC_527/g;s/DC 19/DC_19/g;s/DC 118/DC_118/g;s/SI 144/SI_144/g;s/Dc 120/Dc_120/g;s/Dc 509/Dc_509/g;s/Dc 524/Dc_524/g;s/Dc 122/Dc_122/g;s/Dc 121/Dc_121/g;s/Dc 20/Dc_20/g;s/Dc 112/Dc_112/g;s/Bn 23/Bn_23/g;s/Dc 549/Dc_549/g;s/Bn 32/Bn_32/g;s/DC 532/DC_532/g;s/Bn 35/Bn_35/g;s/Si 21/Si_21/g;s/Dc 551/Dc_551/g;s/Dc 103/Dc_103/g;s/Dc 550/Dc_550/g;s/Dc 105/Dc_105/g;s/Bn 27/Bn_27/g;s/Bn 24/Bn_24/g;s/Bn 174/Bn_174/g;s/Si 16/Si_16/g;s/Dc 500/Dc_500/g;s/Bn 153/Bn_153/g;s/Dc 130/Dc_130/g;s/Bn 629/Bn_629/g;s/Dc 109/Dc_109/g;s/Dc 29/Dc_29/g;s/Bn 185/Bn_185/g;s/Bn 152/Bn_152/g;s/Si 157/Si_157/g;s/Si 136/Si_136/g;s/Si 158/Si_158/g;s/Dc 18/Dc_18/g;s/Si 152/Si_152/g;s/Bn 159/Bn_159/g;s/Si 18/Si_18/g;s/Si 139/Si_139/g;s/Bn 168/Bn_168/g;s/Si 146/Si_146/g;s/Bn 31/Bn_31/g;s/Dc 25/Dc_25/g;s/Si 154/Si_154/g;s/Dc 510/Dc_510/g;s/Si 166/Si_166/g;s/Si 133/Si_133/g;s/Si 135/Si_135/g;s/Dc 124/Dc_124/g;s/Si 141/Si_141/g;s/Si 155/Si_155/g;s/SI 138/SI_138/g;s/BN 155/BN_155/g;s/SI 159/SI_159/g;s/SI 13/SI_13/g;s/SI 17/SI_17/g;s/BN 169/BN_169/g;s/BN 151/BN_151/g;s/BN 154/BN_154/g;s/SI 145/SI_145/g;s/DC 546/DC_546/g;s/SI 149/SI_149/g;s/SI 117/SI_117/g;s/DC 114/DC_114/g;s/DC 27/DC_27/g;s/BN 163/BN_163/g;s/SI 143/SI_143/g;s/DC 540/DC_540/g;s/BN 25/BN_25/g;s/SI 14/SI_14/g;s/54G C/54G_C/g;s/79IA C/79IA_C/g;s/47BR 2343/47BR_2343/g;s/301-9470/301-9470/g;s/DC119/DC119/g;s/8  8 Jan/8__8_Jan/g;s/BR 85/BR_85/g;s/UF 3/UF_3/g;s/BR 121/BR_121/g;s/UF 3/UF_3/g;s/79IA C/79IA_C/g' Univ_of_Florida_Mateescu_BOVF250V1_20220512_FinalReport.txt

cp ../../UF250K_Feb_2024/Illumina/Univ_of_Florida_Mateescu_BOVF250V1_20160331_FinalReport.txt .
cp ../../UF250K_Feb_2024/Illumina/Univ_of_Florida_Mateescu_BOVF250V1_20220505_FinalReport.txt .
cp ../../UF250K_Feb_2024/Illumina/Univ_of_Florida_Mateescu_BOVF250V1_20220512_FinalReport.txt .
