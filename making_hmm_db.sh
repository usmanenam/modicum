cat ncbi_targets.txt | xargs -I {} cat ../../../ncbi_hmms/hmm_PGAP.hmm/{}.HMM >> temp_ncbi_021726.hmm
grep "^ACC" ../../../pfam/Pfam-A.hmm | grep -F -f pfam_targets.txt | awk '{print $2}' > pfam_list_w_versions.txt
hmmfetch -f ../../../pfam/Pfam-A.hmm pfam_list_w_versions.txt > temp_pfam.hmm
cat temp_pfam.hmm temp_ncbi.hmm > pfam_ncbi.hmm
hmmpress pfam_ncbi.hmm