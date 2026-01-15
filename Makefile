.PHONY: m1ad_xnat
	
# M1 S#(SID) actuators to surface displacements frequency response SVD
m1siad:
	cargo r -r --features faer -- \
		-i m1-actuators-segment$(SID) -o m1-segment$(SID)-axial-d \
		--svd -f gmt_m1s$(SID)-actuators_frequency_response.pkl \
		log-space -l 1 -u 100 -n 1000

# M1 S1 actuators to surface displacements frequency response SVD
# for a set of FEM natural frequencies given by their indices
m1s1ad_xnat:
	cargo r -r --features faer -- \
		-i m1-actuators-segment1 -o m1-segment1-axial-d \
		-c --svd --uv u -f gmt_m1s1-actuators_3natural_frequency_response.pkl \
		structural -i 154 -i 414 -i 2045

# M1 S2 actuators to surface displacements frequency response SVD
# for a set of FEM natural frequencies given by their indices
m1s2ad_xnat:
	cargo r -r --features faer -- \
		-i m1-actuators-segment2 -o m1-segment2-axial-d \
		-c --svd --uv u -f gmt_m1s2-actuators_3natural_frequency_response.pkl \
		structural -i 156 -i 413 -i 2047

# M1 S7 actuators to surface displacements frequency response SVD
# for a set of FEM natural frequencies given by their indices
m1s7ad_xnat:
	cargo r -r --features faer -- \
		-i m1-actuators-segment7 -o m1-segment7-axial-d \
		-c --svd --uv u -f gmt_m1s7-actuators_3natural_frequency_response.pkl \
		structural -i 172 -i 441 -i 1894

# M1 actuators to surface displacements frequency response SVD
m1ad:
	cargo r -r --features faer -- \
		-i m1-actuators-segment1 -o m1-segment1-axial-d \
		-i m1-actuators-segment2 -o m1-segment2-axial-d \
		-i m1-actuators-segment3 -o m1-segment3-axial-d \
		-i m1-actuators-segment4 -o m1-segment4-axial-d \
		-i m1-actuators-segment5 -o m1-segment5-axial-d \
		-i m1-actuators-segment6 -o m1-segment6-axial-d \
		-i m1-actuators-segment7 -o m1-segment7-axial-d \
		--svd -f gmt_m1-actuators_frequency_response.pkl \
		log-space -l 2 -u 100 -n 2000
		
# M1 actuators to surface displacements frequency response SVD
# for a set of FEM natural frequencies given by their indices
m1ad_xnat:
	cargo r -r --features faer -- \
		-i m1-actuators-segment1 -o m1-segment1-axial-d \
		-i m1-actuators-segment2 -o m1-segment2-axial-d \
		-i m1-actuators-segment3 -o m1-segment3-axial-d \
		-i m1-actuators-segment4 -o m1-segment4-axial-d \
		-i m1-actuators-segment5 -o m1-segment5-axial-d \
		-i m1-actuators-segment6 -o m1-segment6-axial-d \
		-i m1-actuators-segment7 -o m1-segment7-axial-d \
		-c --svd --uv u -f gmt_m1-actuators_xnatural_frequency_response.pkl \
		structural -i 26 -i 46 -i 161 -i 405 -i 968
		

# M1 wind loads (@ segment vertices) to surface displacements frequency response SVD
m1wd:
	cargo r -r --features faer -- \
		-i ossm1-lcl6-f \
		-o m1-segment1-axial-d \
		-o m1-segment2-axial-d \
		-o m1-segment3-axial-d \
		-o m1-segment4-axial-d \
		-o m1-segment5-axial-d \
		-o m1-segment6-axial-d \
		-o m1-segment7-axial-d \
		--svd -f gmt_m1-wind-loads_frequency_response.pkl \
		log-space -l 2 -u 100 -n 2000
