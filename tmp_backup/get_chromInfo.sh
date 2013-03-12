mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from chromInfo;" sacCer3 | cut -f1-2
