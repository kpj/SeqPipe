########
# Script to set up individual containers for each alignment process
#

set -e
set -u

cur_wd="$(dirname "$0")"
. "$cur_wd/config.sh"

# handle cmdline-arguments
inp_read=
inp_genome_file=
output_dir="$default_output_dir"

usage() {
    echo "Usage: $0 -i <path to read (directory)> -g <path to reference genome file> [-o <output directory>]"
    exit 1
}

while getopts ":i:g:o:" opt; do
    case "${opt}" in
        i)
            inp_read=${OPTARG}
            ;;
        g)
            inp_genome_file=${OPTARG}
            ;;
        o)
            output_dir=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${inp_read}" ] || [ -z "${inp_genome_file}" ] || [ -z "${output_dir}" ]; then
    usage
fi

# process reads
mkdir -p "$output_dir"

for file in $(find "$inp_read" \( -name "*.fastq.gz" -o -name "*.fastq" \)); do
    # book-keeping
    fname=$(basename $file)
    echo -ne "${INFO}" ">> Handling \"$fname\" <<" "${RESET}\n"

    analysis_wd="$output_dir/runs/$fname"
    mkdir -p "$analysis_wd/input"
    cp "$file" "$analysis_wd/input"

    # prepare directory
    if [[ "$fname" == *".gz" ]]; then
        gunzip -c "$analysis_wd/input/$fname" > "$analysis_wd/input/reads.fastq"
    else
        cp "$analysis_wd/input/$fname" "$analysis_wd/input/reads.fastq"
    fi
    cp -r "$cur_wd/scripts/" "$analysis_wd/"
    cp "$inp_genome_file" "$analysis_wd/input/$genome_file"

    # start analysis
    "$cur_wd/pipeline.sh" "$analysis_wd"
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
