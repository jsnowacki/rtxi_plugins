MODULES=IHCs_model

.PHONY: $(MODULES)
all: $(MODULES)

$(MODULES): 
	@$(MAKE) -C $@ 

.PHONY: clean 
clean:
	-for d in $(MODULES); do $(MAKE) -C$$d clean; done
