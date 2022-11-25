#!usr/bin/bash

# convert output name with the correct date

python3 ./sjaracne.py local -e ./input.exp -g ./sig.txt -o ./output_sig_sjaracne_project_2022-08-25_out_.final
python3 ./sjaracne.py local -e ./input.exp -g ./tf.txt -o ./output_tf_sjaracne_project_2022-08-25_out_.final
