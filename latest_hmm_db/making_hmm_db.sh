grep "^ACC" ../../../ncbi_hmms/hmm_PGAP_combined.hmm | grep -F -f ncbi_targets.txt | awk '{print $2}' > ncbi_list_w_versions.txt
hmmfetch -f ../../../ncbi_hmms/hmm_PGAP_combined.hmm ncbi_list_w_versions.txt > temp_ncbi.hmm
grep "^ACC" ../../../pfam/Pfam-A.hmm | grep -F -f pfam_targets.txt | awk '{print $2}' > pfam_list_w_versions.txt
hmmfetch -f ../../../pfam/Pfam-A.hmm pfam_list_w_versions.txt > temp_pfam.hmm
cat temp_pfam.hmm temp_ncbi.hmm > pfam_ncbi.hmm
#may have to change HMM names if the names aren't unique
hmmpress pfam_ncbi.hmm