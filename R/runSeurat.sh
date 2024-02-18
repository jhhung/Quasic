#/bin/bash

iter=$1
input_dir=$2$1
output_dir=$2/R_output
cluster_resolution=$3

echo "---------------------------"
echo "Run Seurat at $1 iteration"
echo "---------------------------"

if [ !  -d $output_dir ]; then
	mkdir $output_dir
fi

# Check requried file
if  [ ! -f "$input_dir/quants_mat.gz" ]; then
	echo "$input_dir/quants_mat.gz is required!"
	echo "ERROR: Failed to run Seurat"
	exit 0
fi
if  [ ! -f "$input_dir/quants_mat_rows.txt" ]; then
	echo "$input_dir/quants_mat_rows.txt is required"
	echo "ERROR: Failed to run Seurat"
	exit 0
fi
if  [ ! -f "$2/quants_mat_cols.txt" ]; then
	echo "$2/quants_mat_cols.txt is required"
	echo "ERROR: Failed to run Seurat"
	exit 0
fi

Rscript R/runSeurat.R $input_dir/quants_mat.gz $iter $output_dir $cluster_resolution

if [ $? -eq 1 ]; then
	echo "ERROR: Something wrong in Seurat"
	exit 0
else
	echo "--------------------------"
	echo "    Run Seurat success    "
	echo "--------------------------"
	exit 1
fi

