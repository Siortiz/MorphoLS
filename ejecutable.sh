python detection.py
chmod 777 sex_seg.sh
./sex_seg.sh
python psf_mask.py
python inputs_galfitm.py
chmod 777 inputs/galfit_20.input
chmod 777 galfitm-1.4.4-linux-x86_64
./galfitm-1.4.4-linux-x86_64 inputs/galfit_20.input
