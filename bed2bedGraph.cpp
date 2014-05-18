// Substitute for bedtools genomecov -i input.bed -g chrom.sizes -bg.
// Runs about 10x faster.
//
// Usage: bed2bedGraph < input.bed > output.bed
//
// Both input and output are zero-based with half-open intervals:
//   chrXYZ start end		- represents (start, end] 
//
// Tomas Kazmar, 2014

#include <cstdio>
#include <cctype>
#include <cstring>
#include <vector>
#include <algorithm>
#include <functional>
#include <limits>

#define INT_LEN (sizeof(int)*8)
#define ITOA_BUF_LEN INT_LEN
#define MAX_PUTS_CACHE_LEN (100000*3*(ITOA_BUF_LEN+1))
#define MAX_CHR_LEN 100
#define MAX_LINE_LEN (MAX_CHR_LEN + 1 + INT_LEN + 1 + INT_LEN + 2000)

char puts_cache[MAX_PUTS_CACHE_LEN+2];
char *puts_cache_pos = puts_cache;
void parse_line(char *line, char *chr, int &start, int &end);
int put_level(char *chr, int from, int to, int level);
void clean_heap(std::vector<int> &heap, char *chr, int &level, int &pos,
	const int start = std::numeric_limits<int>::max());

int main(int argc, char** argv)
{
	char chr[MAX_CHR_LEN+1];		 /// Current chromosome.
	char old_chr[MAX_CHR_LEN+1] = "\0";	 /// Old chromosome.
	std::vector<int> heap;			 /// Heap to store read ends.
	int start;
	int end;

	int level = 0;
	int pos = -1;
	// XXX: read in batches - further speed up?
	char line[MAX_LINE_LEN + 2];
	while (fgets(line, MAX_LINE_LEN, stdin)) {
		parse_line(line, chr, start, end);
		//printf("[Current pos: %d\n", pos);
		
		// Check for new chromosome.
		if (strncmp(chr, old_chr, MAX_CHR_LEN) != 0) {
			// Empty the heap with ends from the old chromosome.
			clean_heap(heap, old_chr, level, pos);
			// Start new chromosome.
			level = 0;
			pos = start;
			strncpy(old_chr, chr, MAX_CHR_LEN);
		} else if (start < pos) {
			fprintf(stderr, "Error: The input file is not sorted. "
				"Line:\n%s\n", line);
			exit(-1);
		}
		// Store the end of the current read.
		heap.push_back(end);
		std::push_heap(heap.begin(), heap.end(), std::greater<int>());

		// Output ranges that end before the start - nothing can change there.
		clean_heap(heap, chr, level, pos, start);

		// Move on and increase level due to current read. (Except when
		// current read starts where the one on the top of the heap ends.
		// In such a case merge the two ranges).
		if (!heap.empty() && heap.front() == start) {
			std::pop_heap(heap.begin(), heap.end(), std::greater<int>());
			heap.pop_back();
		} else {
			if (pos < start) {
				if (level > 0) { // XXX: remove this check to get all zero-valued ranges in between the reads
					pos = put_level(chr, pos, start, level);
				} else {
					pos = start;
				}
			}
			++level;
		}
	}
	// Empty the heap with ends.
	clean_heap(heap, chr, level, pos);
	// Empty the output cache.
	fwrite(puts_cache, 1, puts_cache_pos - puts_cache, stdout);

	return 0;
}


char *itoa_pos_rev(unsigned int val, char* rbuf)
{
    const unsigned int radix = 10;
    char* p = rbuf;
    unsigned int u = val;
    static char digits[] = "0123456789";

    do {
        unsigned int a = u % radix;
        u /= radix;
        *p-- = digits[a];
    } while (u > 0);

    return p;
}


void caching_puts(const char *buf, const size_t len)
{
	if (puts_cache_pos + len > puts_cache + MAX_PUTS_CACHE_LEN) {
		fwrite(puts_cache, 1, puts_cache_pos - puts_cache, stdout);
		puts_cache_pos = puts_cache;
	}
	memcpy(puts_cache_pos, buf, len);
	puts_cache_pos += len;
	*puts_cache_pos++ = '\n';
}


int put_level(char *chr, int from, int to, int level)
{
	static char itoa_buffer[ITOA_BUF_LEN];
	char *itoa_pos;
	itoa_pos = itoa_pos_rev(level, itoa_buffer + ITOA_BUF_LEN);
	*itoa_pos-- = '\t';
	itoa_pos = itoa_pos_rev(to, itoa_pos);
	*itoa_pos-- = '\t';
	itoa_pos = itoa_pos_rev(from, itoa_pos);
	*itoa_pos = '\t';

	const size_t len = strnlen(chr, MAX_CHR_LEN);
	itoa_pos -= len;
	memcpy(itoa_pos, chr, len);
	caching_puts(itoa_pos, itoa_buffer + ITOA_BUF_LEN - itoa_pos + 1);

	return to;
}


void clean_heap(std::vector<int> &heap, char *chr, int &level, int &pos, const int start)
{
	//printf("[start clean_heap(heap, %s, %d, %d, %d)\n", chr, level, pos, start);
	//printf("[ head.front() = %d\n", heap.empty() ? -1 : heap.front());
	while (!heap.empty() && (heap.front() < start)) {
		const int to = heap.front();
		if (pos < to) {
			//printf("[ put_level from clean_heap:\n");
			pos = put_level(chr, pos, to, level);
		}
		--level;
		std::pop_heap(heap.begin(), heap.end(), std::greater<int>());
		heap.pop_back();
	}
	//printf("[finish clean_heap(heap, %s, %d, %d, %d)\n", chr, level, pos, start);
}


int parse_int(char *s, int &val)
{
	char *p = s;
	val = 0;
	while (isdigit(*p)) {
		val *= 10;
		val += *p - '0';
		++p;
	}
	return p - s + 1;
}


void parse_line(char *line, char *chr, int &start, int &end)
{
	const char *delimiters = " \t";
	const char *digits = "0123456789";

	char *pos = line;
	size_t chr_len = strcspn(pos, delimiters);
	if (chr_len == 0) {
		fprintf(stderr, "Error: Zero-length chromosome encountered "
			"while parsing line:\n%s\n", line);
		exit(-1);
	}
	if (chr_len > MAX_CHR_LEN) {
		fprintf(stderr, "Error: Maximum allowed chromosome length exceed "
			"(limit: %d, actual: %d) while parsing line:\n%s\n",
			MAX_CHR_LEN, chr_len, line);
		exit(-1);
	}
	memcpy(chr, pos, chr_len);
	chr[chr_len] = 0;
	pos += chr_len;
	pos = strpbrk(pos, digits);
	if (pos == NULL) {
		fprintf(stderr, "Error parsing line (start): %s", line);
		exit(-1);
	}
	pos += parse_int(pos, start);
	if (!isdigit(*pos)) {
		fprintf(stderr, "Error: Non-numeric characters (%c) in read start. "
			"Line:\n%s\n", *pos, line);
		exit(-1);
	}
	pos = strpbrk(pos, digits);
	if (pos == NULL) {
		fprintf(stderr, "Error: Cannot parse read end. Line: %s", line);
		exit(-1);
	}
	parse_int(pos, end);
	if (!isdigit(*pos)) {
		fprintf(stderr, "Error: Non-numeric characters in read end. "
			"Line:\n%s\n", line);
		exit(-1);
	}
	if (start >= end) {
		fprintf(stderr, "Error: Reads must have a positive length. "
			"Line:\n%s\n", line);
		exit(-1);
	}
	if (start < 0) {
		fprintf(stderr, "Error: Negative read start encountered "
			"while parsing line:\n%s\n", line);
		exit(-1);
	}
	if (end < 0) {
		fprintf(stderr, "Error: Negative read end encountered "
			"while parsing line:\n%s\n", line);
		exit(-1);
	}
}
