#ifndef JABSWORKER_H
#define JABSWORKER_H

#include <QObject>
#include <QStringList>

extern "C" {
#include "script_session.h"
}

/* Runs JaBS script lines in a separate thread (see MainWindow), so that the GUI stays responsive during long simulations and fits.
 * The session must not be accessed by other threads while a run is in progress. */
class JabsWorker : public QObject
{
    Q_OBJECT
public:
    explicit JabsWorker(script_session *session, QObject *parent = nullptr);
    static int runLine(script_session *session, const QString &line); /* Executes one line of script, returns script_command_status */

public slots:
    void run(const QStringList &lines); /* Emits finished() when done */

signals:
    void statusMessage(const QString &msg);
    void finished(bool error, qsizetype lineno);

private:
    script_session *session;
};

#endif // JABSWORKER_H
