#!/bin/bash
# version 2.1.0

function _check_ld_variable {
	local VAR_NAME=${1}
	local VAR_FLAG=${2}
	local LD_PATH=${3}
	
	if [[ -z ${LD_PATH} ]];then
		echo -e "ERROR: user supplied empty value to -${VAR_FLAG}"
		exit 1
	fi
	if [[ $(( $(find ${LD_PATH}[1-9].l2.ldscore.gz -type f 2> /dev/null | wc -l) + $(find ${LD_PATH}1[0-9].l2.ldscore.gz -type f 2> /dev/null | wc -l) + $(find ${LD_PATH}2[0-2].l2.ldscore.gz -type f 2> /dev/null | wc -l) )) -ne 22 ]];then
		echo -e "ERROR: could not find \"${LD_PATH}\"[1-22].l2.ldscore.gz"
		exit 1
	fi
	eval "${VAR_NAME}=${LD_PATH}"
}

function _is_variable_empty {
	local VAR_NAME=${1}
	local VAR_FLAG=${2}
	local VAR_VAL=${3}
	
	if [[ -z ${VAR_VAL} ]];then
		echo -e "ERROR: use -${VAR_FLAG} [option] to supply the ${VAR_NAME}"
		exit 1
	fi
}

function _rename_existing_file {
	local FILE=${1}
	local TIMESTAMP=$(date "+%Y%m%d")
	
	if [[ $(find . -maxdepth 1 -name "${FILE}" -type f 2> /dev/null | wc -l) -ne 0 ]];then
		echo -e "WARNING: directory already contained a file called \"${FILE}\"\n  The existing file was renamed to \"OLD.${TIMESTAMP}.${FILE}\""
		mv ${FILE} OLD.${TIMESTAMP}.${FILE}
	fi
}

function _remove_existing_file {
	local FILE=${1}
	
	if [[ $(find . -maxdepth 1 -name "${FILE}" -type f 2> /dev/null | wc -l) -ne 0 ]];then
		echo -e "WARNING: \"${FILE}\" already existed and was removed"
		rm ./${FILE}
	fi
}

function _construct_matrix {
	local MATRIX_TYPE=${1}
	local CONSTRAIN_VALUE=${2}
	local COLUMN_NUMBER=${3}
	
	for SUMSTAT_NUM in $(seq 1 ${NUM_SUMSTATS});do
		if [[ ${SUMSTAT_NUM} == ${NUM_SUMSTATS} && ${CONSTRAIN_DIAGONAL} == 'T' ]];then
			echo -e "$(tail -n1 sumstats.list).sumstats.gz\t${CONSTRAIN_VALUE}" > new_column
		else
			awk -v COLUMN_NUMBER=${COLUMN_NUMBER} 'BEGIN{OFS="\t"};{print $2,$COLUMN_NUMBER}' < ${SUMSTAT_NUM}.tab > new_column
			if [[ ${CONSTRAIN_DIAGONAL} == 'T' ]];then
				echo -e "$(head -n${SUMSTAT_NUM} < sumstats.list | tail -n1).sumstats.gz\t${CONSTRAIN_VALUE}" | cat - new_column > dum && mv dum new_column
			fi
		fi
		
		if [[ ${SUMSTAT_NUM} == 1 ]];then
			mv new_column ${MATRIX_TYPE}.tmp1
		else
			awk 'BEGIN{OFS="\t"}
			FNR==NR{a[$1]=$2;next}
			FNR!=NR{if($1 in a){print $0,a[$1]}else{print $0}}
			' new_column ${MATRIX_TYPE}.tmp1 > ${MATRIX_TYPE}.tmp2 && mv ${MATRIX_TYPE}.tmp2 ${MATRIX_TYPE}.tmp1
			rm new_column
		fi
	done
	sed -i 's![ \t]*$!!g' ${MATRIX_TYPE}.tmp1
	
	cut -f1 < ${MATRIX_TYPE}.tmp1 > sumstats
	cut -f2- < ${MATRIX_TYPE}.tmp1 > lower_triangle
	for COLUMN in $(seq 1 ${NUM_SUMSTATS});do
		cut -f${COLUMN} < lower_triangle | tail -n+$(( COLUMN + 1 )) > single_column
		paste -sd $'\t' < single_column > single_row
		head -n${COLUMN} lower_triangle | tail -n1 | paste - single_row -d $'\t' >> ${MATRIX_TYPE}.tmp2
		rm single_column single_row
	done
	
	paste sumstats ${MATRIX_TYPE}.tmp2 -d $'\t' | sed 's![ \t]*$!!g' > ${MATRIX_TYPE}.tab
	rm ${MATRIX_TYPE}.tmp? sumstats lower_triangle
}

printf '*%.0s' {1..96}
echo -e '
*
* Script to create a matrix of genetic correlations from a set of GWAS summary statistics
* (C) 2019 Hill F. Ip
* Department of Biological Psychology, Vrije Universiteit Amsterdam / Netherlands Twin Register
* GNU General Public License v3
*'
printf '*%.0s' {1..96}
echo -e '\n'

BEGINTIME=$(date "+%s.%N")
re='^[0-9]+$'
CONSTRAIN_DIAGONAL='F'
SAVE_LOGS='F'
PREFIX=''
INDIR='./'
LDSC=''
REFLD=''
WLD=''
ncors=5
CALL_FLAGS=''

while getopts 'csl:r:w:p:n:i:' flag;do
	case "${flag}" in
		c)	CALL_FLAGS=$(echo "${CALL_FLAGS}; ;c")
			CONSTRAIN_DIAGONAL='T' ;;
		s)	CALL_FLAGS=$(echo "${CALL_FLAGS}; ;s")
			SAVE_LOGS='T' ;;
		n)	CALL_FLAGS=$(echo "${CALL_FLAGS}; ;n")
			PREFIX="${OPTARG}" ;;
		i)	if [[ -z ${OPTARG} ]];then
				echo -e "WARNING: user supplied an empty value to -i\n  Using default value (current directory)"
			elif [[ ${OPTARG: -1} != '/' ]];then
				echo -en "WARNING: value supplied to -i \"${OPTARG}\" did not end on a forward slash\n  Checking whether value is a proper directory..."
				if [[ $(find "${OPTARG}/" -maxdepth 0 -type d 2> /dev/null | wc -l) == 1 ]];then
					echo "done"
					INDIR="${OPTARG}/"
				else
					echo -e "WARNING: \"${OPTARG}\" did not appear to be an existing directory\n  Using default value (current directory)"
				fi
				echo 
			else
				INDIR=${OPTARG}
			fi ;;
		l)	if [[ -z ${OPTARG} ]];then
				echo "ERROR: user supplied empty value to -l"
				exit 1
			fi
			if [[ ${OPTARG: -1} == '/' ]];then
				echo -en "WARNING: the value supplied to -l \"${OPTARG}\" appears to be a directory\n  Looking for \"ldsc.py\" in the directory..."
				if [[ $(find ${OPTARG} -maxdepth 1 -name "ldsc.py" -type f 2> /dev/null | wc -l) == 1 ]];then
					echo "done"
					LDSC="${OPTARG}ldsc.py"
				else
					echo "ERROR: could not find \"ldsc.py\" in \"${OPTARG}\""
					exit 1
				fi
			elif [[ $(echo ${OPTARG} | awk -F'/' '{print $NF}') != "ldsc.py" ]];then
				echo "ERROR: value supplied to -l \"${OPTARG}\" did not point towards \"ldsc.py\""
				exit 1
			else
				LDSC=${OPTARG}
			fi ;;
		r)	_check_ld_variable 'REFLD' 'r' ${OPTARG} ;;
		w)	_check_ld_variable 'WLD' 'w' ${OPTARG} ;;
		p)	if [[ -z ${OPTARG} ]];then
				echo -e "WARNING: user supplied an empty value to -p\n  Using default value (5)"
			elif ! [[ ${OPTARG} =~ ${re} ]];then
				echo -e "WARNING: invalid value supplied to -p \"${OPTARG}\"\n  Using default value (5)"
			elif [[ ${OPTARG} -gt $(nproc) ]];then
				echo -e "WARNING: value supplied to -p (${OPTARG}) is larger than the number of cores on the machine ($(nproc))\n  Using default value (5)"
			else
				ncors=${OPTARG}
			fi ;;
	esac
done

_is_variable_empty 'location of ldsc.py' 'l' ${LDSC}
_is_variable_empty 'reference LD scores' 'r' ${REFLD}
_is_variable_empty 'LD scores to be used as regression weights' 'w' ${WLD}

echo -e "Analysis started at $(date --date @$(echo ${BEGINTIME} | cut -f1 -d'.') "+%a %d-%b-%Y %T %Z")"
echo -en '
Call:
  -l '${LDSC}'
  -r '${REFLD}'
  -w '${WLD}'
  -i '${INDIR}'
  -p '${ncors}
if [[ ! -z ${CALL_FLAGS} ]];then
	echo " \\"
	echo "${CALL_FLAGS}" | sed 's!^; ;!!g' | sed 's!; ;! \n!g' | sed 's!^!  -!g'
	echo ''
else
	echo -e "\n"
fi
for FILE in $(echo "rg.tab rg_se.tab cti.tab cti_se.tab collected_results.tab");do
	_rename_existing_file ${FILE}
done

if [[ $(find . -maxdepth 1 -name "*.rg.log" -type f 2> /dev/null | wc -l) -ne 0 ]];then
	echo -e "WARNING: directory already contained files ending on \".rg.log\". Some of these files might have been overwritten during this run"
fi

for FILE in $(echo "sumstats.list ldsc_errors.log rg.list dum.rg.list");do
	_remove_existing_file ${FILE}
done

if [[ ${INDIR} != './' ]];then
	if [[ $(find . -maxdepth 1 -name "*.sumstats.gz" -type f 2> /dev/null | wc -l) -ne 0 ]];then
		echo "WARNING: user supplied a different input-directory, but the current directory contained munged summary statistics. These will be included in the rg-matrix"
	fi
	if [[ $(find ${INDIR} -maxdepth 1 -name "*.sumstats.gz" -type f 2> /dev/null | wc -l) == 0 ]];then
		echo -e "ERROR: the supplied input-directory (\"${INDIR}\") did not contain munged summary statistics\n  These files should end on \".sumstats.gz\""
		exit 1
	fi
	echo -en "Copying sumstats from ${INDIR} ..."
	cp ${INDIR}*.sumstats.gz .
	echo "done"
fi

echo -en "Listing sumstats..."
ls *.sumstats.gz 2> /dev/null | sed 's!.sumstats.gz!!g' > sumstats.list
NUM_SUMSTATS=$(wc -l < sumstats.list)
if [[ ${NUM_SUMSTATS} == 0 ]];then
	echo -e "ERROR: could not find any munged summary statistics\n  These files should end on \".sumstats.gz\""
	exit 1
elif [[ ${NUM_SUMSTATS} == 1 ]];then
	echo -e "ERROR: directory contained exactly one munged summary statistics file"
	exit 1
fi

for SUMSTAT_NUM in $(seq 1 ${NUM_SUMSTATS});do
	if [[ ${CONSTRAIN_DIAGONAL} == 'T' ]];then
		if [[ ${SUMSTAT_NUM} == ${NUM_SUMSTATS} ]];then
			break
		else
			FIRST_CORRELATION_SUMSTAT_NUM=$(( SUMSTAT_NUM + 1 ))
		fi
	else
		FIRST_CORRELATION_SUMSTAT_NUM=${SUMSTAT_NUM}
	fi
	CURRENT_SUMSTAT=$(awk -v SUMSTAT_NUM=${SUMSTAT_NUM} 'NR==SUMSTAT_NUM' < sumstats.list)
	for CORRELATION_SUMSTAT_NUM in $(seq ${FIRST_CORRELATION_SUMSTAT_NUM} ${NUM_SUMSTATS});do
		echo -e "${SUMSTAT_NUM}\t${CORRELATION_SUMSTAT_NUM}\t${CURRENT_SUMSTAT}" >> dum.rg.list
	done
	tail -n+${FIRST_CORRELATION_SUMSTAT_NUM} sumstats.list | paste dum.rg.list - >> rg.list
	rm dum.rg.list
done
echo -e "done\n  ${NUM_SUMSTATS} sumstats were found"

echo -en "Running pairwise LDSC (NOTE: this might take a while)..."
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
echo "done"

if [[ $(ls tmp_*.rg.log 2> /dev/null | wc -l) -ne $(wc -l < rg.list) ]];then
	echo -e "ERROR: $(ls tmp_*.rg.log 2> /dev/null | wc -l)/$(wc -l < rg.list) LDSC analyses ran"
	exit 1
fi

echo -en "Collecting LDSC results..."
echo -e "p1\tp2\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se" > collected_results.tab
while read SUMSTAT_NUM CORRELATION_SUMSTAT_NUM SUMSTAT1 SUMSTAT2;do
	if [[ $(grep -E "^p1 +p2" tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log | wc -l) == 0 ]];then
		echo -e "ERROR: the analysis between ${SUMSTAT1} and ${SUMSTAT2} failed to return estimates of the genetic correlation (other analyses might also have failed)"
		exit 1
	fi
	tail -n+$(( $(grep -En "^p1 +p2" tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log | cut -f1 -d":") + 1 )) < tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log | head -n1 | tee -a ${SUMSTAT_NUM}.tab collected_results.tab > /dev/null
done < rg.list
echo "done"

if [[ $(grep -iE 'warning|error' tmp_*.rg.log | wc -l) -ne 0 ]];then
	echo -e "WARNING: found warning(s)/error(s) in the LDSC-logs. All logs will be saved"
	SAVE_LOGS='T'
fi

if [[ ${SAVE_LOGS} == 'T' ]];then
	while read SUMSTAT_NUM CORRELATION_SUMSTAT_NUM SUMSTAT1 SUMSTAT2;do
		mv tmp_${SUMSTAT_NUM}_${CORRELATION_SUMSTAT_NUM}.rg.log ${SUMSTAT1}_${SUMSTAT2}.rg.log
	done < rg.list
else
	rm tmp_*.rg.log
fi

echo -en "Constructing the matrices..."
_construct_matrix 'rg' 1 3
_construct_matrix 'rg_se' 0 4
_construct_matrix 'cti' 1 11
_construct_matrix 'cti_se' 0 12
echo "done"

if [[ ! -z ${PREFIX} ]];then
	for FILE in $(echo "rg rg_se cti cti_se collected_results");do
		_rename_existing_file ${PREFIX}${FILE}.tab
		mv ${FILE}.tab ${PREFIX}${FILE}.tab
	done
fi
echo -e "Files are saved as: ${PREFIX}rg.tab, ${PREFIX}rg_se.tab, ${PREFIX}cti.tab, ${PREFIX}cti_se.tab, and collected_results.tab"

rm sumstats.list rg.list
for SUMSTAT_NUM in $(seq 1 ${NUM_SUMSTATS});do
	rm ${SUMSTAT_NUM}.tab 2> /dev/null
done

ENDTIME=$(date "+%s.%N")
TOTAL_TIME=$(date -u -d "0 $ENDTIME sec - ${BEGINTIME} sec" "+%H:%M:%S.%3N")
echo -e "\nScript finished at $(date --date @$(echo ${ENDTIME} | cut -f1 -d'.') "+%a %d-%b-%Y %T %Z")\nElapsed time: ${TOTAL_TIME}"
