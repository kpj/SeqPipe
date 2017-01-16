########
# Script to align reads against many genomes in succession
# Genome directory should have the following structure:
#   00-genome01.fa, 10-genome01.fa, 20-genome01.fa, ...
#

set -e
set -u

if [[ $# -ne 2 && $# -ne 3 ]]; then
    echo "Usage: $0 <path to read (directory)> <path to genome directory> [output directory]"
    exit
fi

cur_wd="$(dirname "$0")"
inp_read="$1"
inp_genome_dir="$2"

. "$cur_wd/../config.sh"

output_dir="${3:-sequential_results}"
mkdir -p "$output_dir"

# helper functions
function do_mapping {
    local  __retvar=$4

    local cur_read="$1"
    local genome_file="$2"
    local stage="$3"

    local map_dir="$output_dir/$stage/mapping"
    local out_dir="$output_dir/$stage/output"

    mkdir -p "$map_dir"
    local map_log="$map_dir/log.txt"

    "$cur_wd/../main.sh" "$cur_read" "$genome_file" "$map_dir" &> "$map_log"

    # copy bam files (as fastq) to output directory
    mkdir -p "$out_dir"
    imp_file="not_aligned_reads.bam" # continue only with unaligned reads

    for read_file in $(find "$cur_read" -name "*.fastq"); do
        read_base="$(basename $read_file)"

        path="$map_dir/runs/$read_base/$imp_file"
        out_file="$out_dir/${read_base}" # ${imp_file%%.bam}_

        bedtools bamtofastq -i "$path" -fq "$out_file"

        # additional files
        EXTRA_FILES=("aligned_reads.bam")

        for file in ${EXTRA_FILES[@]}; do
            path="$map_dir/runs/$read_base/$file"
            out_file="$out_dir/${read_base}_extra/${file%%.bam}.fastq"

            mkdir -p "$(dirname $out_file)"
            bedtools bamtofastq -i "$path" -fq "$out_file"
        done
    done

    eval $__retvar="'$out_dir'"
}

# process reads
stage="initial"

echo -ne "${BACKGROUND}"
echo "Preparing..."
mkdir -p "$output_dir/$stage/input/"
for file in $(find "$inp_read" -name "*.fastq.gz"); do
    out="$output_dir/$stage/input/$(basename ${file%%.gz})"
    echo " > $file -> $out"
    gunzip -c "$file" > "$out"
done
echo "Starting..."
echo -ne "${RESET}"

for file in $(find "$inp_genome_dir" -name "*.fa" | sort -h); do
    fname=$(basename $file)
    echo -ne "${INFO}" ">> $fname ($stage) <<" "${RESET}\n"

    next_stage="${stage}_${fname%%.fa}"
    read_input="$output_dir/$stage/input/"

    # store some extra information
    echo "$fname" > "$output_dir/$stage/info.txt"

    #echo "Input: \"$read_input\""
    do_mapping "$read_input" "$file" "$stage" read_output
    #echo "Output: \"$read_output\""

    # put current stage's output as next stage's input
    mkdir -p "$output_dir/$next_stage/input"
    cp "$read_output/"* "$output_dir/$next_stage/input/" &> /dev/null || true

    stage="$next_stage"
done
