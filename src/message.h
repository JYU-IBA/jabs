/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_MESSAGE_H
#define JABS_MESSAGE_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef enum jabs_msg_level {
    MSG_DEBUG = 0,
    MSG_VERBOSE = 1,
    MSG_INFO = 2,
    MSG_IMPORTANT = 3,
    MSG_WARNING = 4,
    MSG_ERROR = 5
} jabs_msg_level;

#define JABS_DEFAULT_VERBOSITY (MSG_INFO)

extern jabs_msg_level jabs_message_verbosity;

void jabs_message_printf(jabs_msg_level level, FILE *f, const char *format, ...);
void jabs_message(jabs_msg_level level, const char *format, ...);
const char *jabs_message_level_str(jabs_msg_level level);
#ifdef __cplusplus
}
#endif
#endif // JABS_MESSAGE_H
