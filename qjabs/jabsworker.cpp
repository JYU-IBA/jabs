#include "jabsworker.h"
extern "C" {
#include "script.h"
#include "script_command.h"
#include "generic.h"
#include "message.h"
#include "defaults.h"
}

JabsWorker::JabsWorker(script_session *session, QObject *parent)
    : QObject(parent), session(session)
{
}

int JabsWorker::runLine(script_session *session, const QString &line)
{
    if(line.isEmpty()) {
        return SCRIPT_COMMAND_SUCCESS;
    }
    jabs_message(MSG_INFO, "%s%s\n", PROMPT, qPrintable(line));
    return script_execute_command(session, qPrintable(line));
}

void JabsWorker::run(const QStringList &lines)
{
    qsizetype lineno = 0;
    bool error = false;
    for(const QString &line : lines) {
        lineno++;
        if(jabs_abort_requested()) {
            jabs_message(MSG_ERROR, "Run aborted by user.\n");
            error = true;
            break;
        }
        QString command_str = line.split(" ").first();
        const script_command *c = script_command_find(session->commands, qPrintable(command_str));
        if(c && c->f == script_fit) {
            emit statusMessage(QString("Running a fit."));
        } else if(c && c->f == script_simulate) {
            emit statusMessage(QString("Running a simulation."));
        }
        if(runLine(session, line) < 0) {
            error = true;
            break;
        }
        if(session->file_depth > 0) {
            if(script_process(session) < 0) {
                error = true;
                break;
            }
        }
    }
    emit finished(error, lineno);
}
