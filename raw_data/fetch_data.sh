# fetch all the sra records
prefetch --option-file accession_list.txt
# convert them into fastqs
cat accession_list.txt | xargs fasterq-dump --outdir "raw_fastqs/"
# remove .sra records
rm -fr SRR*/
# gzip the fastq files
for FILE in raw_fastqs/*.fastq;
do
    gzip ${FILE}
done
