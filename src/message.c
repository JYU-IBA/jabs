#include <stdio.h>
#include <stdarg.h>
#include "message.h"

jabs_msg_level jabs_message_verbosity;
static const char *jabs_msg_levels[MSG_ERROR+1] = {"Debug", "Verbose", "Default", "Important", "Warning", "Error"};

void jabs_message(jabs_msg_level level, const char * restrict format, ...) {
    if(level < jabs_message_verbosity) {
        return;
    }
    va_list argp;
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
}

void jabs_message_printf(jabs_msg_level level, FILE *f, const char * restrict format, ...) {
    if(f == stderr && level < jabs_message_verbosity) {
        return;
    }
    va_list argp;
    va_start(argp, format);
    vfprintf(f, format, argp);
    va_end(argp);
}

const char *jabs_message_level_str(jabs_msg_level level) {
    if(level <= MSG_ERROR) {
        return jabs_msg_levels[level];
    }
    return NULL;
}
