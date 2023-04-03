export

GCC_VER_TO_INT = $(shell echo $(1) | awk -F. '{print $$3+100*($$2+100*$$1)}')
GCC_VERSION := $(shell gcc -dumpfullversion)
#GCC_MAJOR_VERSION := $(shell echo $(GCC_VERSION) | awk -F. '{print $$1}')
GCC_MAJOR_VERSION := $(firstword $(subst ., , $(GCC_VERSION)))
GCC_VERSION_INT := $(call GCC_VER_TO_INT, $(GCC_VERSION))
GCC_VER_GTE10 := $(shell [ $(GCC_VERSION_INT) -ge 100000 ] || echo 0 )

#It will get called during every expansion phase, printing lots of error messages
#GCC_VER_CMP = $(shell test $(GCC_VERSION_INT) -$(1) $(call GCC_VER_TO_INT, $(2) ) || echo 0 )

# Ex:
# ifeq ($(GCC_VER_GTE10),)
# $(info True)
# else
# $(info False)
# endif

