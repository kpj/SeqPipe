# Script archive

At the end of each alignment process, all scripts will be executed on the resulting (sorted) bam file. They should thus have the following signature:
```
<script> <path to bam>
```
The currently supported script formats are:
* `*.py`: Python scripts
* `*.sh`: Bash scripts