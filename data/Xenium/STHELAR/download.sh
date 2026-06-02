FILENAME="sdata_lung_s1.zarr.zip"
DATA_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/${FILENAME%%.*}"
mkdir -p $DATA_PATH
cd $DATA_PATH
curl -L -O "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/146/S-BIAD2146/Files/STHELAR/sdata_slides/${FILENAME}"
unzip $FILENAME