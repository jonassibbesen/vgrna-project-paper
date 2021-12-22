set -e

SCRIPT=${1}
ARGS=${2}
JOB_NAME=${3}
JOB_INFO=${4}

S3_DIR=$(dirname "${SCRIPT}")
SCRIPT_BASE=$(basename "${SCRIPT}")
SCRIPT_PREFIX="${SCRIPT_BASE%.sh}"

LOG_NAME=${JOB_NAME}_log.txt

bash -c "
	echo ---
	echo ${JOB_NAME}
	echo ${JOB_INFO}
	echo ${SCRIPT}
	echo ${ARGS}
	echo ---
	aws s3 cp ${SCRIPT} .
	env ${ARGS} /bin/sh ${SCRIPT_BASE}
	echo ---
" 2>&1 | tee ${LOG_NAME}

aws s3 cp ${LOG_NAME} ${S3_DIR}/
