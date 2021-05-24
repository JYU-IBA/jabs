#include <stdlib.h>
#include <string.h>

char **string_to_argv(const char *str) { /* TODO: move this generic function to elsewhere! */
    char *s = strdup(str);
    char *s_split = s;
    s[strcspn(s, "\r\n")] = 0;
    size_t len = strlen(s);
    size_t n = 0;
    char *col;
    while((col = strsep(&s_split, " \t")) != NULL) {
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
    for(pos = 0; pos < len; pos++) {
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
