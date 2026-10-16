/* Keep checks active in Release builds, including calls that build the fixture. */
#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "freesasa.h"

int main(void)
{
    /* Export three calculations assembled through both public tree APIs. */
    freesasa_structure *s = freesasa_structure_new();
    assert(s);
    assert(freesasa_structure_add_atom(s, " CA ", "ALA", "1", 'A', 0, 0, 0) == FREESASA_SUCCESS);
    freesasa_result *r = freesasa_calc_structure(s, NULL);
    assert(r);
    freesasa_node *tree = freesasa_tree_init(r, s, "first");
    freesasa_node *other = freesasa_tree_init(r, s, "joined");
    assert(tree && other);
    assert(freesasa_tree_add_result(tree, r, s, "added") == FREESASA_SUCCESS);
    assert(freesasa_tree_join(tree, &other) == FREESASA_SUCCESS);
    assert(other == NULL);
    for (int skip = 0; skip < 2; ++skip) {
        FILE *out = tmpfile();
        assert(out);
        assert(freesasa_tree_export(out, tree, FREESASA_RSA | (skip ? FREESASA_OUTPUT_SKIP_REL : 0)) == FREESASA_SUCCESS);
        rewind(out);
        char line[512];
        int headers = 0, residues = 0, totals = 0, chains = 0, names = 0;
        while (fgets(line, sizeof line, out)) {
            headers += strncmp(line, "REM  Absolute and relative", 26) == 0;
            residues += strncmp(line, "RES ", 4) == 0;
            totals += strncmp(line, "TOTAL ", 6) == 0;
            chains += strncmp(line, "CHAIN", 5) == 0;
            names += strstr(line, "SASAs for first") != NULL || strstr(line, "SASAs for added") != NULL || strstr(line, "SASAs for joined") != NULL;
        }
        assert(headers == 3 && residues == 3 && totals == 3 && chains == 3 && names == 3);
        fclose(out);
    }
    /* Linux's always-full device verifies propagation of buffered write errors. */
    FILE *full = fopen("/dev/full", "w");
    if (full) {
        assert(freesasa_tree_export(full, tree, FREESASA_RSA) == FREESASA_FAIL);
        fclose(full);
    }
    freesasa_node_free(tree);
    freesasa_result_free(r);
    freesasa_structure_free(s);
    puts("RSA multi-result, joined-tree, skip-relative and output-error checks passed.");
}
