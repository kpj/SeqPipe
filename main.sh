########
# Script to set up individual containers for each alignment process
#

set -e
set -u


if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <path to read directory> <path to reference genome file>"
    exit
fi

cur_wd="$(dirname "$0")"
inp_read_dir="$1"
inp_genome_file="$2"

. "$cur_wd/config.sh"

# process reads
mkdir -p "$output_dir"

for file in $(find "$inp_read_dir" -name "*.fastq.gz"); do
    # book-keeping
    id=$(basename $file | cut -d'_' -f1)
    fname=$(basename $file)

    analysis_wd="$output_dir/runs/$id"
    mkdir -p "$analysis_wd/input"
    cp "$file" "$analysis_wd/input"

    # prepare directory
    gunzip -c "$analysis_wd/input/$fname" > "$analysis_wd/input/reads.fastq"
    cp -r "$cur_wd/scripts/" "$analysis_wd/"
    cp "$inp_genome_file" "$analysis_wd/input/$genome_file"

    # start analysis
    "$cur_wd/map_reads.sh" "$analysis_wd"
done


# gather results
if [[ ! -d "$output_dir/runs/" ]]; then
    echo "No results were generated..."
    exit -1
fi

res_dir="$output_dir/results"
mkdir -p "$res_dir"

echo "Gathering results"
for dir in "$output_dir/runs/"*; do
    id="$(basename $dir)"
    echo " > $id"
    mkdir -p "$res_dir/$id/"
    for res in "$dir/results/"*; do
        echo "  - $res"
        cp "$res" "$res_dir/$id/"
    done
done
