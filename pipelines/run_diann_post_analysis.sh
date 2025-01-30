#!/bin/bash

help_message () {
	echo ""
	echo "Options:"
	echo ""
	echo "	-o PATH	Analysis output path"
	echo "	-n STRING	Screen name. Used for output zip. Defaults to out.zip"
  echo "	-i PATH	DIANN result path"
  echo "	-e FILE	experimental information TSV file: file<tab>sample<tab>sample_name<tab>condition<tab>replicate"
  echo "	-c STRING	experimental comparisons list: comp1 comp2 comp2"
  echo ""
  echo "  -h --help	Display help"

	echo "";
}

export EXE_PATH=$(dirname "$0")
in=""
basepath="";
name="out"
exp=""
comps=""

# load in params
# load in params
SHORT_OPTS="h:i:o:e:c:n:"
LONG_OPTS='help,in,out,name,exp,comps'

OPTS=$(getopt -o $SHORT_OPTS --long $LONG_OPTS -- "$@")
# make sure the params are entered correctly
if [ $? -ne 0 ];
then
    help_message;
    exit 1;
fi

while true; do
    # echo $1
	case "$1" in
      -h | --help) help_message; exit 1; shift 1;;
      -i) in=$2; shift 2;;
      -o) basepath=$2; shift 2;;
      -n) name=$2; shift 2;;
      -e) exp=$2; shift 2;;
      -c) comps=$2; shift 2;;
      --) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done

if [ "$in" = "" ]; then
    echo "Please provide a diann result path."
    help_message; exit 1
else
    echo "## DIANN result path: $in"
fi

if [ "$basepath" = "" ]; then
    echo "Please provide an output base path"
    help_message; exit 1
else
    echo "## Output basepath: ${basepath}"
    echo "## Screen name: ${name}"
fi

if [ "$exp" = "" ]; then
    echo "Please provide an experimental design TSV file: file<tab>sample<tab>sample_name<tab>condition<tab>replicate"
    help_message; exit 1
else
    echo "## Experimental info TSV: ${basepath}"
fi

if [ "$comps" = "" ]; then
    echo "Please provide an experimental comparisons list (space seperated): comp1 comp2 comp2"
    help_message; exit 1
else
    echo "## Comparisons: ${comps}"
fi

echo "loading modules"
module load StdEnv/2023 r/4.4.0 netcdf/4.9.2 gdal/3.9.1

# filter protein groups to kepp proteotypic protein groups only
echo "generating proteoptypic protein group matrix"
perl -ne '
chomp($_);
my @t = split("\t",$_);
my @prot_ident = split(";",$t[0]);
if(scalar(@prot_ident) == 1){
  print $_ . "\n";
}
' ${in}/report.pg_matrix.tsv > ${in}/report.pg_matrix.proteoptypic.tsv

# run fragpipe
echo "running fragpipe"

Rscript ${EXE_PATH}/../R/fragpipe_analysis.R \
--matrix ${in}/report.pg_matrix.proteoptypic.tsv \
--out ${basepath}/fragpipe \
--design ${exp}

# extract columns based on comparisons
echo "generate fragpipe filtered results based on comparisons to perform"
comps_regex=$(echo "$comps" | sed -r 's/[ ]+/\|/g')
awk -v get="Protein.Group|Genes|First.Protein.Description|$comps_regex" '
BEGIN{
  FS=OFS="\t"
}
FNR==1{
  for(i=1;i<=NF;i++)
    if($i~get)cols[++c]=i
}
{
  for(i=1; i<=c; i++)
    printf "%s%s", $(cols[i]), (i<c ? OFS : ORS)
}' ${basepath}/fragpipe/report.pg_matrix.proteoptypic.dge.tsv > ${basepath}/fragpipe/report.pg_matrix.proteoptypic.dge.filtered.tsv

# gen reactome
echo "Reactome analysis"
for comp in $comps
do
  for geneset in KEGG Reactome MSigDB
  do
      echo "############## running ${comp} using geneset ${geneset} and PIN ${pin}"
      Rscript ${EXE_PATH}/../R/gen_reactome.R \
      -i ${basepath}/fragpipe/report.pg_matrix.proteoptypic.dge.tsv \
      -o ${basepath}/reactome \
      -s mouse -g Genes \
      -c ${comp} \
      --gene_source ${geneset}
  done
done

# gen gsea
echo "GSEA analysis"
for comp in $comps
do
  for geneset in kegg go
  do
    echo "running ${comp} using ${geneset} gene set"
    Rscript ${EXE_PATH}/../R/gen_gsea.R \
    -m ${basepath}/fragpipe/report.pg_matrix.proteoptypic.dge.tsv \
    -c ${comp} -s Mm -g ${geneset} \
    -o ${basepath}/gsea
  done
done

echo "Wrapup and zip results"
cd ${basepath}
mkdir -p diann_reports
cp ${in}/*matrix.tsv diann_reports/
cp ${in}/report.stats.tsv diann_reports/
cp ${in}/report.log.txt diann_reports/
cp ${in}/report-lib.tsv.speclib diann_reports/
cp ${in}/report-fasta-database.fasta diann_reports/

t=$(basename ${exp})
if [ ! -f ${basepath}/${t} ]; then
    cp ${exp} .
fi
zip -r ${name}.zip diann_reports fragpipe gsea reactome ${t}
cd -