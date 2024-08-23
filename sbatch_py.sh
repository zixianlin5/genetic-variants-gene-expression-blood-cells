#!/bin/bash
#SBATCH --job-name=paga_analysis       # Job name
#SBATCH --output=paga_analysis-%j.out  # Standard output file
#SBATCH --error=paga_analysis-%j.err   # Standard error file
#SBATCH --mem=200G                     # Memory allocation
#SBATCH --mail-type=END,FAIL           # Send email at end of job or on job failure
#SBATCH --mail-user=email@example.com   # Your email address

hostname

LOCAL_TMP_DIR=/tmp/$SLURM_JOB_ID
mkdir -p $LOCAL_TMP_DIR

LOCAL_OUTPUT_DIR=$LOCAL_TMP_DIR/pyoutput/$SLURM_JOB_ID
mkdir -p $LOCAL_OUTPUT_DIR
echo "Output will be saved to: $LOCAL_OUTPUT_DIR"

cp -r /path/to/data/pyinput $LOCAL_TMP_DIR/
cp -r /path/to/data/scanpy_env $LOCAL_TMP_DIR/
cp -r /path/to/data/py_script $LOCAL_TMP_DIR/

cd $LOCAL_TMP_DIR

source $LOCAL_TMP_DIR/scanpy_env/bin/activate

python $LOCAL_TMP_DIR/py_script/pyscript_name.py --local_tmp_dir $LOCAL_TMP_DIR --local_output_dir $LOCAL_OUTPUT_DIR

# Copy the output folder back to the shared file system
if [ "$(ls -A $LOCAL_OUTPUT_DIR)" ]; then
    echo "The following files will be copied from $LOCAL_OUTPUT_DIR to /path/to/data/pyoutput/:"
    
    for file in $LOCAL_OUTPUT_DIR/*; do
        echo "$(basename $file)"
    done

    cp -r $LOCAL_OUTPUT_DIR /path/to/data/pyoutput/
    
    echo "Files have been copied to /path/to/data/pyoutput/"
else
    echo "No files to copy from $LOCAL_OUTPUT_DIR"
fi

rm -rf $LOCAL_TMP_DIR