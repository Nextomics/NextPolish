TOP_DIR := $(PWD)
BIN_DIR := $(TOP_DIR)/bin
SUB_DIR := lib util

TARGER_UTIL := seq_split seq_count bwa samtools minimap2

define cp_b
if [ -e $(TOP_DIR)/util/$(N)/$(N) ]; then \
	cp $(TOP_DIR)/util/$(N)/$(N) $(BIN_DIR); \
elif [ -e $(TOP_DIR)/util/$(N) ]; then \
	cp $(TOP_DIR)/util/$(N) $(BIN_DIR); \
else \
	cp $(TOP_DIR)/lib/$(N) $(BIN_DIR); \
fi;
endef

all: make_bin
	$(foreach N, $(SUB_DIR), make -C $(N);)
	@$(foreach N, $(TARGER_UTIL), $(call cp_b, $(N)))

make_bin:
ifeq ($(wildcard $(BIN_DIR)),)
	mkdir $(BIN_DIR)
endif

clean:
	rm -rf $(BIN_DIR)
	$(foreach N, $(SUB_DIR), make -C $(N) clean;)
	
