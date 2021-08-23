#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "win_compat.h"

char **string_to_argv(const char *str) {
    char *s = strdup(str);
    char *s_split = s;
    s[strcspn(s, "\r\n")] = 0;
    size_t len = strlen(s);
    size_t n = 0;
    char *col;
    while((col = strsep(&s_split, " \t")) != NULL) { /* TODO: parse quotation marks so that 'foo "bar baz"' becomes out[0] == "foo", out[1] == "bar baz"; */
        if(*col == '\0') {
            continue;
        }
        n++;
    }
    if(!n) {
        free(s);
    }
    char **out = malloc(sizeof(char *) * (n + 1));
    out[0] = s;
    size_t pos;
    size_t i = 1;
    for(pos = 0; pos < len && i < n; pos++) {
        if(s[pos] == '\0' /*&& s[pos-1] != '\0'*/) {
            out[i] = s+pos+1;
            if(*out[i] != '\0') { /* Consecutive delimeters (turned to '\0' by strsep above) are ignored */
                i++;
            }
        }
    }
    out[n] = NULL;
    return out;
}

char *argv_to_string(int argc, char * const *argv) {
    if(argc < 1)
        return NULL;
    char *s = NULL;
    size_t len = 0;
    for(int i = 0; i < argc; i++) {
        len += strlen(argv[i]);
        len++;
    }
    s = malloc(sizeof(char) * len);
    *s = '\0';
    char *sp = s;
    for(int i = 0; i < argc; i++) {
        strcat(sp, argv[i]);
        sp += strlen(argv[i]);
        *sp = ' ';
        sp++;
    }
    *sp = '\0';
    return s;
}

FILE *fopen_file_or_stream(const char *filename, const char *mode) {
    FILE *f;
    if(!filename) {
        if(*mode == 'r')
            f = stdin;
        else
            f = stderr;
    } else if(strlen(filename) == 1 && *filename == '-')
        f = stdout;
    else {
        f = fopen(filename, mode);
    }
    if(!f && filename) {
        fprintf(stderr, "Can not open file \"%s\" (mode %s)\n", filename, mode);
    }
    return f;
}

void fclose_file_or_stream(FILE *f) {
    if(f == stdout || f == stderr || f == stdin)
        return;
    fclose(f);
}

char *strdup_non_null(const char *s) {
    if(s)
        return strdup(s);
    else
        return NULL;
}
