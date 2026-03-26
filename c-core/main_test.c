#include <stdio.h>
#include "dna_utils.h"
#include "needleman_wunsch.h"
#include "suffix_tree.h"
#include "trie_errors.h"
#include "segment_tree.h"
#include "interval_tree.h"
#include "union_find.h"

int main() {
    printf("GenomeX C-Core Test Suite\n");
    printf("==========================\n");
    printf("This is a stub test runner over the integration of multiple structures.\n");
    printf("Once all .c files are implemented, you should expand this test.\n");
    printf("Testing dna_utils: validation -> %d\n", dna_validate("ATGC"));
    return 0;
}
