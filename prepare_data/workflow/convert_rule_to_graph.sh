SNAKEFILE=$1
RULE=$2
UNTIL=${3:-$RULE}
snakemake --dag -s $SNAKEFILE -U $UNTIL -n $RULE | sed 1d | dot -Tsvg > graph.svg