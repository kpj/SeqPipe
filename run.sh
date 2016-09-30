########
# Script to set up individual containers for each alignment process
#

set -e
set -u


if [[ $# -ne 2 ]] ; then
    echo "Usage: $0 <path to reads> <path to reference genome file>"
    exit
fi

rel_dir="$(dirname "$0")"
read_dir="$1"
genome_file="$2"

# process reads
wd="mapping_result"
mkdir -p "$wd"

for file in $(find "$read_dir" -name "*.fastq.gz"); do
    # book-keeping
    id=$(basename $file | cut -d'_' -f1)
    fname=$(basename $file)

    cwd="$wd/runs/$id"
    mkdir -p "$cwd"
    cp "$file" "$cwd"

    # prepare directory
    gunzip -c "$cwd/$fname" > "$cwd/reads.fastq"
    cp "$rel_dir/scripts/"* "$cwd"
    cp "$genome_file" "$cwd"

    # start analysis
    "$rel_dir/map_reads.sh" "$cwd"
done


# gather results
img_dir="$wd/images"
mkdir -p "$img_dir"

for dir in "$wd/runs/"*; do
    for img in "$dir/images/"*; do
        echo $img
        cp "$img" "$img_dir"
    done
done
