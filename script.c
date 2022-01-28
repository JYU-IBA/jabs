/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <jibal_units.h>
#include "fit.h"
#include "generic.h"
#include "spectrum.h"
#include "script.h"
#include "jabs.h"
#include "message.h"

int script_process(script_session *s, const char *filename) {
    if(script_session_load_script(s, filename)) {
        return EXIT_FAILURE;
    }
    static const char *prompt = PROMPT;
    script_command_status status = SCRIPT_COMMAND_SUCCESS;
    while(s->file_depth > 0) {
        script_file *sfile = s->files[s->file_depth - 1];
        int interactive = sfile->interactive;
        if(status == SCRIPT_COMMAND_EXIT || (status != SCRIPT_COMMAND_SUCCESS && !interactive)) { /* on exit, close all nested script files. When non-interactive, close scripts on error until interactive (or exit). */
            script_file_close(sfile);
            s->file_depth--;
            continue;
        }
        if(interactive) {
            fputs(prompt, stderr);
        }
        if(script_file_getline(sfile) > 0) {
            if(!interactive) {
                jabs_message(MSG_INFO, stderr, "%s%s\n", prompt, sfile->line);
            }
            status = script_execute_command(s, sfile->line);
            if(!interactive && status != SCRIPT_COMMAND_SUCCESS) {
                jabs_message(MSG_ERROR, stderr, "Error %i (%s) on line %zu in file \"%s\". Aborting.\n", status,
                             script_command_status_to_string(status), sfile->lineno, sfile->filename);
            }
        } else {
            status = SCRIPT_COMMAND_EOF;
            if(interactive)
                status = SCRIPT_COMMAND_EXIT;
        }
    }
    if(s->files[0]->interactive) {
        jabs_message(MSG_INFO, stderr, "\nBye!\n");
    }
    fflush(stderr);
    return status;
}
