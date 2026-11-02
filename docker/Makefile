IMAGE ?= brentp/genome-flow-align:0.8.0

.PHONY: align
align:
	docker build -t $(IMAGE) -f align.Dockerfile ..

.PHONY: glnexus
glnexus:
	docker build -t brentp/glnexus:1.4.5 -f glnexus.Dockerfile ..
