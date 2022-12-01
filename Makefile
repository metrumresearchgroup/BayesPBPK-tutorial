build-image:
	docker build . -t bayesianpbpk:latest --no-cache

# Run with  docker run -it --rm  -p 8443:8443 -p 8787:8787 bayespbpk:latest
