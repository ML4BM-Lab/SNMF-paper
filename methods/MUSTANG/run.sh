#!/bin/bash
# run.sh
# Usage: bash run.sh /path/to/data

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    echo "Usage: bash submit_pipeline.sh /path/to/data output/path"
    exit 1
fi

OUTPUT_PATH=$2
if [ -z "$OUTPUT_PATH" ]; then
    echo "Usage: bash submit_pipeline.sh /path/to/data output/path"
    exit 1
fi

TECHNOLOGY=$3
if [ -z "$TECHNOLOGY" ]; then
    echo "Usage: bash submit_pipeline.sh /path/to/data output/path"
    exit 1
fi

# Step 1
jid1=$(sbatch --parsable ./load_data/load_data.sbatch $DATA_PATH $OUTPUT_PATH)

# Step 2
jid2=$(sbatch --parsable --dependency=afterok:$jid ./clustering/clustering.sbatch $OUTPUT_PATH)

# Step 3
jid3=$(sbatch --parsable --dependency=afterok:$jid2 ./graph/graph.sbatch $DATA_PATH $OUTPUT_PATH)

# Step 4
jid4=$(sbatch --parsable --dependency=afterok:$jid3 ./deconvolution/deconvolution.sbatch $DATA_PATH $OUTPUT_PATH $TECHNOLOGY)
