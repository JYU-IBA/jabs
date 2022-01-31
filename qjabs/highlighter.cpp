/* Original copyright notice of highlighter.cpp below, used
 * under the terms of the BSD license.
 *
 * Modifications by Jaakko Julin for JaBS. Modified version may be distributed
 * under the terms of the GPL (version 2 or later) provided that the terms of
 * the BSD license are also met. Note that the modified version may not be
 * distributed under the BSD license or Qt commercial license. */

/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "highlighter.h"
extern "C" {
#include <../generic.h>
}

Highlighter::Highlighter(QTextDocument *parent)
    : QSyntaxHighlighter(parent)
{
    singleLineCommentFormat.setForeground(Qt::gray);
    commandFormat.setForeground(Qt::darkBlue);
    commandFormat.setFontWeight(QFont::Bold);
    variableFormat.setForeground(Qt::darkYellow);
    variableFormat.setFontWeight(QFont::StyleItalic);
}

void Highlighter::setSession(const script_session *session)
{
    Highlighter::session = session;
}

void Highlighter::highlightBlock(const QString &text)
{
    qsizetype len = text.length();
    if(len == 0) {
        setCurrentBlockState(0);
        return;
    }
    if(text.at(0) == '#') {
        setFormat(0, text.length(), singleLineCommentFormat);
        setCurrentBlockState(0);
        return;
    }
    int argc = 0;
    char **argv = string_to_argv(qPrintable(text), &argc);
    if(!argv) {
        return;

    }
    highlightArgv(argc, argv);
    argv_free(argv, argc);


    return;

    for (const HighlightingRule &rule : qAsConst(highlightingRules)) {
        QRegularExpressionMatchIterator matchIterator = rule.pattern.globalMatch(text);
        while (matchIterator.hasNext()) {
            QRegularExpressionMatch match = matchIterator.next();
            setFormat(match.capturedStart(), match.capturedLength(), rule.format);
        }
    }
    setCurrentBlockState(0);
}

void Highlighter::highlightArgv(int argc, char **argv) {
    const script_command *cmds = session->commands;
    const script_command *c_parent = NULL;
    const char *argv_start = argv[0];
    while(argc && cmds) {
        const script_command *c = script_command_find(cmds, argv[0]);
        if(!c) {
            return;
        }
        while(c) { /* Subcommand found */
            size_t arg_len = strlen(argv[0]);
            if(c->f) {
                setFormat(argv[0] - argv_start, arg_len, commandFormat);
            } else if(c->var) {
                setFormat(argv[0] - argv_start, arg_len, variableFormat);
            }
            if(c->subcommands) {
                setFormat(argv[0] - argv_start, arg_len, commandFormat);
                argc--;
                argv++;
                cmds = c->subcommands;
                c_parent = c;
                break;
            } else {
                argc--;
                argv++;
                c = NULL;
            }
        }
    }
}
