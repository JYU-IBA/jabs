#include <stdio.h>
#include <stdarg.h>
#include "message.h"

void jabs_message(jabs_msg_level level, FILE *f, const char * restrict format, ...) {
    va_list argp;
    va_start(argp, format);
    vfprintf(f, format, argp);
    va_end(argp);
}
