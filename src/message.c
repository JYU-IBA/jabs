#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
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

static int jabs_progress_len = 0; /* Length of progress line currently shown (0 if none) */

void jabs_progress(double fraction, const char *format, ...) { /* Command line: progress line (with fraction as percentage) is updated in place, only if stderr is a terminal */
    if(jabs_message_verbosity > MSG_INFO || !isatty(fileno(stderr))) {
        return;
    }
    va_list argp;
    va_start(argp, format);
    fputc('\r', stderr);
    int len = 0;
    if(fraction >= 0.0) {
        len += fprintf(stderr, "[%5.1lf%%] ", 100.0 * (fraction > 1.0 ? 1.0 : fraction));
    }
    len += vfprintf(stderr, format, argp);
    va_end(argp);
    if(len < jabs_progress_len) { /* Clear the rest of previous (longer) line */
        fprintf(stderr, "%*s", jabs_progress_len - len, "");
    }
    jabs_progress_len = len > 0 ? len : 0;
    fflush(stderr);
}

void jabs_progress_done(void) {
    if(jabs_progress_len > 0) {
        fprintf(stderr, "\r%*s\r", jabs_progress_len, "");
        fflush(stderr);
        jabs_progress_len = 0;
    }
}
