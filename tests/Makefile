TESTS=starts_where_ends\
 overlap\
 one_bp_read\
 one_bp_overlap\
 fail_empty_read\
 fail_too_long_chromosome\
 fail_unsorted\
 fail_incomplete_read\
 fail_negative_length\
 fail_bad_start_pos_read\
 fail_negative_start

.PHONY: all

B2BG=../bed2bedGraph

all: $(TESTS)

starts_where_ends:
	./expect-success.sh $@ "$(B2BG)" $@.bed $@.bg

overlap:
	./expect-success.sh $@ "$(B2BG)" $@.bed $@.bg

one_bp_read:
	./expect-success.sh $@ "$(B2BG)" $@.bed $@.bg

one_bp_overlap:
	./expect-success.sh $@ "$(B2BG)" $@.bed $@.bg

fail_empty_read:
	./expect-failure.sh $@ "$(B2BG)" $@.bed $@.bg

fail_too_long_chromosome:
	./expect-failure.sh $@ "$(B2BG)" $@.bed

fail_incomplete_read:
	./expect-failure.sh $@ "$(B2BG)" $@.bed

fail_negative_length:
	./expect-failure.sh $@ "$(B2BG)" $@.bed

fail_bad_start_pos_read:
	./expect-failure.sh $@ "$(B2BG)" $@.bed

fail_negative_start:
	./expect-failure.sh $@ "$(B2BG)" $@.bed

fail_unsorted:
	./expect-failure.sh $@ "$(B2BG)" $@.bed
