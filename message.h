/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef _JABS_MESSAGE_H_
#define _JABS_MESSAGE_H_

#include <stdio.h>

typedef enum jabs_msg_level {
    MSG_DEBUG = 0,
    MSG_VERBOSE = 1,
    MSG_INFO = 2,
    MSG_WARNING = 3,
    MSG_ERROR = 4
} jabs_msg_level;

void jabs_message(jabs_msg_level level, FILE *f, const char *format, ...);
#endif // _JABS_MESSAGE_H_
