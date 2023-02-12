#include <stdio.h>
#include <stdarg.h>
#include "message.h"

jabs_msg_level jabs_message_verbosity;

void jabs_message(jabs_msg_level level, FILE *f, const char * restrict format, ...) {
    if(level < jabs_message_verbosity) {
        return;
    }
    va_list argp;
    va_start(argp, format);
    vfprintf(f, format, argp);
    va_end(argp);
}
