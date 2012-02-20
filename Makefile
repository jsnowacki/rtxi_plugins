# Here list all the plugin subfolders separated by a space
MODULES=IHCs_model IHCs_model_ICs HH_neuron_ICs 

# This construct allows a usage of make with multiple job option (-j#)
.PHONY: $(MODULES)
all: $(MODULES)

$(MODULES): 
	@$(MAKE) -C $@ 

.PHONY: clean 
clean:
	-for d in $(MODULES); do $(MAKE) -C$$d clean; done

.PHONY: install 
install:
	-for d in $(MODULES); do $(MAKE) -C$$d install; done
