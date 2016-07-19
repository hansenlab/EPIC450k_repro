# Esteller data:

cd /dcs01/hansen/hansen_lab1/meth_epic_esteller
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75073/suppl/GSE75073_Signals_450k.txt.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75073/suppl/GSE75073_Signals_Epic.txt.gz

# Illumina data:
cd /dcs01/hansen/hansen_lab1/meth_epic_illumina
curl -O ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-demo-dataset.zip
curl -O http://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-sample-sheet.zip


# ENCODE 450k data:
cd /dcs01/hansen/hansen_lab1/meth_450k_encode
curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/supplemental/wgEncodeHaibMethyl450BetaValues.txt