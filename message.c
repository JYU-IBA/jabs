#include <stdio.h>
#include <stdarg.h>
#include "message.h"

void jabs_message(jabs_msg_level level, const char * restrict format, ...) {
    va_list argp;
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
}
