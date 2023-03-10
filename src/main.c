/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See LICENSE.txt for the full license.

 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <jibal.h>
#include <jibal_cs.h>

#include "jabs_debug.h"
#include "defaults.h"
#include "generic.h"
#include "git.h"
#include "options.h"
#include "sample.h"
#include "simulation.h"
#include "spectrum.h"
#include "fit.h"
#include "script.h"
#include "script_session.h"
#include "script_command.h"
#include "idf2jbs.h"
#include "message.h"

int idf2jbs(int argc, char * const *argv) {
    char *filename_out = NULL;
    if(argc == 1) {
        idf_error idferr = idf_parse(argv[0], &filename_out);
        if(idferr == IDF2JBS_SUCCESS) {
            jabs_message(MSG_INFO, "Success. Wrote script to file \"%s\"\n", filename_out);
            free(filename_out);
            return EXIT_SUCCESS;
        } else {
            jabs_message(MSG_ERROR, "IDF2JBS failed with error code %i (%s).\n", idferr, idf_error_code_to_str(idferr));
            return EXIT_FAILURE;
        }
    } else {
        jabs_message(MSG_INFO, "Usage (idf2jbs): jabs idf2jbs <idf file>\nNote that idf file must have suffix .xml or .idf.\nOn successful run .jbs (simulation script) and .dat (spectrum) files will be created.\n");
        return EXIT_FAILURE;
    }
}

int main(int argc, char * const *argv) {
    DEBUGMSG("This is a debug build, version %s", jabs_version());
    if(git_populated()) {
        DEBUGMSG("%sGit build. branch %s, commit %s, date %s",
                 git_dirty() ? "Dirty " : "", git_branch(), git_commit_sha1(), git_commit_date());
    }
    if(argc >= 2 && strcmp(argv[1], "idf2jbs") == 0) {
        argc -= 2;
        argv += 2;
        return idf2jbs(argc, argv);
    }
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        jabs_message(MSG_ERROR, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return EXIT_FAILURE;
    }
    script_session *session = script_session_init(jibal, NULL);
    if(!session) {
        jabs_message(MSG_ERROR, "Can not initialize session.\n");
        jibal_free(jibal);
        return EXIT_FAILURE;
    }
    cmdline_options *cmd_opt = cmdline_options_init();
    read_options(cmd_opt, &argc, &argv);
    jabs_message_verbosity = cmd_opt->verbose;
    DEBUGMSG("Verbosity %i", jabs_message_verbosity);
    if(argc == 0) {
        cmd_opt->interactive = TRUE;
    }
    if(!cmd_opt->interactive) {
        greeting(FALSE);
    }
    int status = 0;
    for(int i = 0; i < argc; i++) {
        if(script_session_load_script(session, argv[i])) {
            return EXIT_FAILURE;
        }
        status = script_process(session);
        DEBUGMSG("Script %i/%i given from command line has been processed. Status: %s", i + 1, argc, script_command_status_to_string(status));
        if(status != SCRIPT_COMMAND_SUCCESS) {
            return status;
        }
    }
    DEBUGSTR("All scripts given from command line have been processed.");
    if(cmd_opt->interactive) {
        greeting(TRUE);
        if(script_session_load_script(session, NULL)) {
            return EXIT_FAILURE;
        }
        status = script_process(session);
    }
    script_session_free(session);
    cmdline_options_free(cmd_opt);
    jibal_free(jibal);
    return status == SCRIPT_COMMAND_FAILURE ? EXIT_FAILURE : EXIT_SUCCESS;
}
