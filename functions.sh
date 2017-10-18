#!/usr/bin/env bash

set -eu


function aggregate_samtools_stats {
    local bam_files=( "$@" )

    for fname in "${bam_files[@]}"; do
        echo "> $fname"
        samtools flagstat "$fname"
        echo
    done
}
