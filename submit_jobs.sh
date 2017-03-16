#!/bin/bash 
 source ~/.bashrc

# training RNAs
# rnas="1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 1R7W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 2Y95 "

# testing RNAs
rnas="1SCL 1Z2J 2L3E 2M24 2MHI 2N6S 2NBZ 1ZC5 2LUN 2M4W 2MNC 2N6T 2NC0 2M5U 2MXL 2N6W 2NCI 1UUU 2JWV 2M12 2M8K 2N2O 2N6X 2QH2 5A17 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 5A18 2M22 2MFD 2N4L 2NBY 5KQE 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 1R7W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 2Y95"
# loop over rnas and get data
for rna in ${rnas}
do
	echo "sed 's/XXXX/${rna}/g' scripts/job.sh | sed 's/YYYY/larmord/g'" | bash > scripts/job_tmp.sh
	sbatch scripts/job_tmp.sh

	echo "sed 's/XXXX/${rna}/g' scripts/job.sh | sed 's/YYYY/ramsey/g'" | bash > scripts/job_tmp.sh
	sbatch scripts/job_tmp.sh
done



