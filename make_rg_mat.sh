#!/bin/bash

#####   README   #####
##
## This script calculates pairwise genetic correlations between a set of GWAS summary statistics (version 2.0.0)
## Last edited: 01-MAR-2019
##
## Input: "munged" GWAS summary statistics, see https://github.com/bulik/ldsc
## NOTE: script further requires Python version >=2.7 & <3; preferably 2.7.9
##
## Output:
##  3x tab-delmited text files; N rows by (N+1) columns:
##   - rg.tab:	symmetric matrix of genetic correlations
##   - se.tab:	symmetric matrix of SE's for the genetic correlations
##   - cti.tab:	symmetric matrix of cross-trait-intercepts
##  1x tab-delimited text file containing all estimated values reported by LDSC; (N+1) rows (including header) by 12 columns
##
## How to use:
##  1) copy this script and all sumstats to the working directory
##  2) run the script using the following command: bash script.make_rg_mat.sh
## NOTE: the script assumes no files ending on ".sumstats.gz" are present except for the GWASs for which the genetic correlations are to be calculated
##
## Additional flags:
##  -c			-- restrict the diagonal of both rg.tab and cti.tab to all 1's, and restrict the diagonal of se.tab to 0's
##  -s			-- save all LDSC log-files
##  -n [value]	-- add prefix to output files
##  -p [value]	-- change the number of parallel runs of LDSC to start; default: 5
##  -l [value]	-- supply the location of ldsc.py (see https://github.com/bulik/ldsc ); default: /home/vubiopsy/GWAS/resources/software/ldsc/ldsc.py (NOTE: access to this file might be restricted)
##  -r [value]	-- indicate which files to use as reference LDscores (see https://github.com/bulik/ldsc ); default: /home/vubiopsy/GWAS/resources/LDScores/eur_w_ld_chr/ (NOTE: access to this file might be restricted)
##  -w [value]	-- indicate which files to use as regression weights (see https://github.com/bulik/ldsc ); default: /home/vubiopsy/GWAS/resources/LDScores/eur_w_ld_chr/ (NOTE: access to this file might be restricted)
##
## Example command to run the script with additional flags:
##  bash script.make_rg_mat.sh -scp 8 -l ~/ldsc/ldsc.py
##
#####   end README   #####

#####   Prepare workspace for script   #####

## define functions

## function to check whether the variable points towards EXACTLY 22 files ending on "[1-22].l2.ldscore.gz"
## function requires two arguments (in this order):
##  (1) the location of the LD scores (e.g. "~/LDScores/baseline." --> function looks for "~/LDScores/baseline.[1-22].l2.ldscore.gz")
##  (2) the name of the variable (i.e. "REFLD" or "WLD")
function _check_ld_variable {
	if [[ $(( $(find "${1}"[1-9].l2.ldscore.gz 2> /dev/null | wc -l) + $(find "${1}"1[0-9].l2.ldscore.gz 2> /dev/null | wc -l) + $(find "${1}"2[0-2].l2.ldscore.gz 2> /dev/null | wc -l) )) -ne 22 ]];then
		echo -e 'WARNING: could not find '${1}'[1-22].l2.ldscore.gz\n  Returning to default: /home/vubiopsy/GWAS/resources/LDScores/eur_w_ld_chr/ (NOTE: access to this directory might be restricted)'
	else
		eval "${2}=${1}"
	fi
}

## function receives a file name as argument and checks whether the file exists in the current directory. If so, a warning-message is displayed and the file is removed
function _remove_existing_file {
	if [[ $(find . -maxdepth 1 -name ${1} 2> /dev/null | wc -l) -ne 0 ]];then
		echo -e 'WARNING: found "'${1}'". This was removed'
		rm ./${1}
	fi
}

## function to create a full symmetric matrix
## function takes three arguments (in this order):
##  (1) type of data to be put in the matrix (also used as name);
##  (2) the default value (for CONSTRAIN_DIAGONAL=T); and
##  (3) in which column the data are found
## Input: tab-delimited file without header, with the following name structure: <number>.tab (e.g. 1.tab, 2.tab, etc.)
function _construct_matrix {
	## define variables
	MATRIX_TYPE=${1}
	CONSTRAIN_VALUE=${2}
	COLUMN_NUMBER=${3}

	## construct lower triangle of the matrix
	for SUMSTAT_NUM in $(seq 1 ${TOT_NUM_SUMSTATS});do
		if [[ ${SUMSTAT_NUM} == ${TOT_NUM_SUMSTATS} && ${CONSTRAIN_DIAGONAL} == 'T' ]];then
			echo -e "$(tail -n1 sumstats.list)\t${CONSTRAIN_VALUE}" > new_column
		else
			awk -v COLUMN_NUMBER=${COLUMN_NUMBER} 'BEGIN{OFS="\t"};{print $2,$COLUMN_NUMBER > "new_column"}' < ${SUMSTAT_NUM}.tab
			if [[ ${CONSTRAIN_DIAGONAL} == 'T' ]];then
				echo -e "$(head -n${SUMSTAT_NUM} < sumstats.list | tail -n1)\t${CONSTRAIN_VALUE}" | cat - new_column > dum && mv dum new_column
			fi
		fi

		if [[ ${SUMSTAT_NUM} == 1 ]];then
			mv new_column ${MATRIX_TYPE}.tmp1
		else
			join -a1 -t $'\t' -1 1 -2 1 ${MATRIX_TYPE}.tmp1 -o auto new_column > dum && mv dum ${MATRIX_TYPE}.tmp1
			rm new_column
		fi
	done
	sed -i 's![ \t]*$!!g' ${MATRIX_TYPE}.tmp1

	## transpose lower triangle to get the upper triangle and combine the two to make the full symmetric matrix
	cut -f1 < ${MATRIX_TYPE}.tmp1 > sumstats
	cut -f2- < ${MATRIX_TYPE}.tmp1 > lower_triangle
	for COLUMN in $(seq 1 ${TOT_NUM_SUMSTATS});do
		cut -f${COLUMN} < lower_triangle | tail -n+$((COLUMN + 1)) > single_column
		paste -sd $'\t' < single_column > single_row
		head -n${COLUMN} lower_triangle | tail -n1 | paste - single_row -d $'\t' >> ${MATRIX_TYPE}.tmp2
		rm single_column single_row
	done
	paste sumstats ${MATRIX_TYPE}.tmp2 -d $'\t' | sed 's![ \t]*$!!g' > ${MATRIX_TYPE}.tab
	rm ${MATRIX_TYPE}.tmp? sumstats lower_triangle
}

## initialize variables (including default values)
BEGINTIME=$(date +%s.%N)
re='^[0-9]+$'
CONSTRAIN_DIAGONAL=''
SAVE_LOGS=''
PREFIX=''
LDSC='/home/vubiopsy/GWAS/resources/software/ldsc/ldsc.py'
REFLD='/home/vubiopsy/GWAS/resources/LDScores/eur_w_ld_chr/'
WLD='/home/vubiopsy/GWAS/resources/LDScores/eur_w_ld_chr/'
ncors=5
CALL_FLAGS=''

## if flags are set, update corresponding variables (after some checks)
while getopts 'csl:r:w:p:n:' flag;do
	case "${flag}" in
		c)	CALL_FLAGS=$(echo "${CALL_FLAGS}; ;c")
			CONSTRAIN_DIAGONAL='T' ;;
		s)	CALL_FLAGS=$(echo "${CALL_FLAGS}; ;s")
			SAVE_LOGS='T' ;;
		n)	CALL_FLAGS=$(echo "${CALL_FLAGS}; ;n ${OPTARG}")
			PREFIX=${OPTARG} ;;
		l)	if [[ $(echo ${OPTARG:(-1)}) == '/' ]];then
				echo -en 'WARNING: the value supplied to -l "'${LDSC}'" appears to be a directory\n  Looking for "ldsc.py" in the directory...'
				if [[ $(find ${OPTARG} -maxdepth 1 -name ldsc.py 2> /dev/null | wc -l) == 1 ]];then
					echo 'done'
					LDSC=${OPTARG}ldsc.py
				else
					echo -e 'WARNING: could not find "ldsc.py"\n  Returning to default: /home/vubiopsy/GWAS/resources/software/ldsc/ldsc.py (NOTE: access to this file might be restricted)'
				fi
			elif [[ $(echo ${OPTARG} | awk -F'/' '{print $NF}') != 'ldsc.py' ]];then
				echo -e 'WARNING: the value supplied to -l ('${LDSC}') does not end with "ldsc.py"\n  Returning to default: /home/vubiopsy/GWAS/resources/software/ldsc/ldsc.py (NOTE: access to this file might be restricted)'
			else
				LDSC=${OPTARG}
			fi ;;
		r)	_check_ld_variable ${OPTARG} 'REFLD' ;;
		w)	_check_ld_variable ${OPTARG} 'WLD' ;;
		p)	if ! [[ ${OPTARG} =~ ${re} ]];then
				echo -e 'WARNING: invalid value supplied to -p '${ncors}'\n  Using default value of 5'
			elif [[ ${OPTARG} -gt $(nproc) ]];then
				echo -e 'WARNING: value supplied to -p ('${ncors}') is larger than the number of cores in the machine ('$(nproc)')\n  Using default value of 5'
			else
				## supplied value is both a number, and less than or equal to the number of cores of the machine/node
				ncors=${OPTARG}
			fi ;;
	esac
done

## print opening text
printf '*%.0s' {1..96}
echo -e '
*
* Script to create a matrix of genetic correlations from a set of GWAS summary statistics'
echo -e '* (C) 2019 Hill F. Ip
* Department of Biological Psychology, Vrije Universiteit Amsterdam / Netherlands Twin Register
* GNU General Public License v3
*'
printf '*%.0s' {1..96}
echo -e "\nCall:\n  -l ${LDSC} \\"
echo -e "  -r ${REFLD} \\"
echo -e "  -w ${WLD} \\"
echo -en "  -p ${ncors}"
if [[ ! -z ${CALL_FLAGS} ]];then
	echo " \\"
	echo "${CALL_FLAGS}" | sed 's!^; ;!!g' | sed 's!; ;! \\\n!g' | sed 's!^!  -!g'
	echo ''
else
	echo -e '\n'
fi
echo -e 'Analysis started at '$(date --date @$(echo ${BEGINTIME} | cut -f1 -d'.') +"%a %d-%b-%Y %T %Z")

## remove/rename files from previous runs
for FILE in $(echo "rg.tab se.tab cti.tab collected_results.tab");do
	if [[ $(find . -maxdepth 1 -name ${FILE} 2> /dev/null | wc -l) -ne 0 ]];then
		echo -e 'WARNING: directory already contains a file called "'${FILE}'"\n  The existing file was renamed to: OLD.'$(date +%Y%m%d)'.'${FILE}''
		mv ${FILE} OLD.$(date +%Y%m%d).${FILE}
	fi
done

if [[ $(find *.rg.log -maxdepth 1 2> /dev/null | wc -l) -ne 0 ]];then
	echo -e 'WARNING: directory contains files ending on ".rg.log" and might have been overwritten during the current run'
fi

_remove_existing_file 'sumstats.list'
_remove_existing_file 'ldsc_errors.log'
_remove_existing_file 'rg.list'

#####   Run bivariate LDSC   #####

echo -en 'Listing sumstats...'
ls *.sumstats.gz 2> /dev/null | sed 's!.sumstats.gz!!g' > sumstats.list
TOT_NUM_SUMSTATS=$(wc -l < sumstats.list)
if [[ ${TOT_NUM_SUMSTATS} == 0 ]];then
	echo -e 'ERROR: directory contains no input files'
	exit 1
elif [[ ${TOT_NUM_SUMSTATS} == 0 ]];then
	echo -e 'ERROR: directory contains only one input file'
	exit 1
fi

## make list detailing which sumstats need to be correlated with which (since only values for the lower triangle of the matrices need to be calculated)
for SUMSTAT_NUM in $(seq 1 ${TOT_NUM_SUMSTATS});do
	if [[ ${CONSTRAIN_DIAGONAL} == 'T' ]];then
		if [[ ${SUMSTAT_NUM} == ${TOT_NUM_SUMSTATS} ]];then
			break
		else
			FIRST_CORRELATION_SUMSTAT_NUM=$((SUMSTAT_NUM + 1))
		fi
	else
		FIRST_CORRELATION_SUMSTAT_NUM=${SUMSTAT_NUM}
	fi
	CURRENT_SUMSTAT=$(awk -v SUMSTAT_NUM=${SUMSTAT_NUM} 'NR==SUMSTAT_NUM' < sumstats.list)
	for CORRELATION_SUMSTAT_NUM in $(seq ${FIRST_CORRELATION_SUMSTAT_NUM} ${TOT_NUM_SUMSTATS});do
		echo -e "${SUMSTAT_NUM}\t${CORRELATION_SUMSTAT_NUM}\t${CURRENT_SUMSTAT}" >> dum.rg.list
	done
	tail -n+${FIRST_CORRELATION_SUMSTAT_NUM} sumstats.list | paste dum.rg.list - >> rg.list
	rm dum.rg.list
done

echo -e "done\n  ${TOT_NUM_SUMSTATS} files were found"

## run LDSC
echo -en 'Running pairwise LDSC (NOTE: this might take a while)...'
while read SUMSTAT_NUM CORRELATION_SUMSTAT_NUM SUMSTAT1 SUMSTAT2;do
	((k=k%ncors)); ((k++==0)) && wait
	{
		python ${LDSC} \
		--rg ${SUMSTAT1}.sumstats.gz,${SUMSTAT2}.sumstats.gz \
		--ref-ld-chr ${REFLD} \
		--w-ld-chr ${WLD} \
		--out tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg > /dev/null 2>&1
	} &
done < rg.list
wait
echo 'done'

if [[ $(ls tmp_*.rg.log 2> /dev/null | wc -l) -ne $(wc -l < rg.list) ]];then
	echo -en 'ERROR: '
	if [[ $(ls tmp_*.rg.log 2> /dev/null | wc -l) == 0 ]];then
		echo -en '0'
	else
		echo -en 'only '$(ls tmp_*.rg.log 2> /dev/null | wc -l)
	fi
	echo -e '/'$(wc -l < rg.list)' LDSC analyses ran'
	exit 1
fi

## collect estimates into tables; 
echo -en 'Collecting LDSC results...'
echo -e 'p1\tp2\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se' > collected_results.tab
while read SUMSTAT_NUM CORRELATION_SUMSTAT_NUM SUMSTAT1 SUMSTAT2;do
	if [[ $(grep -E '^p1 +p2'< tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log | wc -l) == 0 ]];then
		echo 'ERROR: the analysis between '${SUMSTAT1}' and '${SUMSTAT2}' failed to return estimates of the genetic correlation (other analyses might also have failed)'
		exit 1
	fi
	tail -n+$(( $(grep -En '^p1 +p2' < tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log | cut -f1 -d':') + 1 )) < tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log | head -n1 | tee -a ${SUMSTAT_NUM}.tab collected_results.tab >/dev/null
done < rg.list
echo 'done'

## look for warning(s)/error(s) in the LDSC logs
grep -iEn 'warning|error' tmp_*.rg.log > ldsc_errors.log
if [[ $(wc -l < ldsc_errors.log) -ne 0 ]];then
	echo -e '  WARNING: found warning(s)/error(s) in the LDSC-logs. All logs will be saved'
	SAVE_LOGS='T'
else
	rm ldsc_errors.log
fi
## check if LDSC logs need to be saved; if so, rename the logs, if not, remove them
if [[ ${SAVE_LOGS} == 'T' ]];then
	while read SUMSTAT_NUM CORRELATION_SUMSTAT_NUM SUMSTAT1 SUMSTAT2;do
		mv tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log ${SUMSTAT1}_${SUMSTAT2}.rg.log
	done < rg.list
else
	rm tmp_*.rg.log
fi

#####   Construct matrices   #####

echo -en 'Constructing the matrices...'

_construct_matrix 'rg' 1 3
_construct_matrix 'se' 0 4
_construct_matrix 'cti' 1 11

echo 'done'

#####   Finish script   #####

## if requested, rename output-files
if [[ ! -z ${PREFIX} ]];then
	for FILE in $(echo "rg.tab se.tab cti.tab collected_results.tab");do
		if [[ $(find . -maxdepth 1 -name ${PREFIX}${FILE} 2> /dev/null | wc -l) -ne 0 ]];then
			echo -e 'WARNING: directory already contains a file called "'${PREFIX}${FILE}'"\n  The existing file was renamed to: OLD.'$(date +%Y%m%d)'.'${PREFIX}${FILE}
			mv ${PREFIX}${FILE} OLD.$(date +%Y%m%d).${PREFIX}${FILE}
		fi
		mv ${FILE} ${PREFIX}${FILE}
	done
fi
echo -e "Files are saved as: ${PREFIX}rg.tab, ${PREFIX}se.tab, ${PREFIX}cti.tab and ${PREFIX}collected_results.tab\n"

## clear workspace
rm sumstats.list rg.list
for SUMSTAT_NUM in $(seq 1 ${TOT_NUM_SUMSTATS});do
	rm ${SUMSTAT_NUM}.tab 2> /dev/null
done

## print closing text
ENDTIME=$(date +%s.%N)
TOTAL_TIME=$(date -u -d "0 $ENDTIME sec - $BEGINTIME sec" +"%H:%M:%S.%3N")
echo -e 'Analysis finished at '$(date --date @$(echo ${ENDTIME} | cut -f1 -d'.') +"%a %d-%b-%Y %T %Z")
echo "Elapsed time: ${TOTAL_TIME}"
