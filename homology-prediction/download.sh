# Download the Ensembl106 compara gene trees and the corrected SCORPiOs trees

DOWNLOAD_DIR=data/download

mkdir -p $DOWNLOAD_DIR
cd $DOWNLOAD_DIR

curl -O http://ftp.ensembl.org/pub/misc/aqua-faang/ensembl_trees/gene_trees/Compara.106.protein_default.nhx.emf.gz 
curl -O http://ftp.ensembl.org/pub/misc/aqua-faang/ensembl_trees/scorpio_corrected_trees/106_PostSCORPiOs.tgz
# Note: We ended up not using the SCORPiOS corrected trees because it sometimes
# wrongly tries to force LORe gene-trees into AORe topology.
# However we need to download it because we use the species tree that is included.

tar xvf 106_PostSCORPiOs.tgz
