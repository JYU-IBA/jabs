#ifndef COLORPROFILE_H
#define COLORPROFILE_H

#include <QColor>
#include <QIcon>
#include <QMap>
#include <QSettings>
#include <QString>
#include <QStringList>
#include <QVector>

/* Color profiles of the GUI. Two built-in profiles (light and dark) and user defined profiles (stored in QSettings).
 * The selected profile can be "Automatic", which chooses the light or dark profile based on the color scheme of the OS.
 * An invalid color means "system default" (palette of the OS), which is allowed for some roles (e.g. backgrounds). */

struct ColorRole {
    QString key;
    QString label;
    QString group;
    QColor light; /* Default in light profile (invalid = system) */
    QColor dark; /* Default in dark profile (invalid = system) */
    bool systemAllowed; /* Can be set to system default */
};

class ColorProfile
{
public:
    static const QString automaticName;
    static const QString lightName;
    static const QString darkName;
    static const int n_reaction_colors = 10;

    ColorProfile() = default;
    ColorProfile(const QString &name, bool builtin);
    QString name;
    bool builtin = false;
    QColor color(const QString &key) const; /* Invalid if system default */
    QColor color(const QString &key, const QColor &systemColor) const; /* systemColor if system default */
    void setColor(const QString &key, const QColor &color);
    QColor reactionColor(int index) const;
    Qt::ColorScheme colorScheme() const; /* Color scheme (light/dark) of the application for this profile, Unknown means follow the OS */

    static const QVector<ColorRole> &roles();
    static ColorProfile light();
    static ColorProfile dark();
    static ColorProfile automatic();
    static bool systemIsDark();
    static bool isBuiltinName(const QString &name); /* Includes "Automatic" */
    static QStringList customProfileNames(QSettings &settings);
    static ColorProfile load(QSettings &settings, const QString &name); /* Built-in profiles or custom ones. Unknown name gives automatic. */
    static ColorProfile active(QSettings &settings); /* Profile selected in settings, "Automatic" resolved */
    static void saveCustomProfiles(QSettings &settings, const QList<ColorProfile> &profiles); /* Replaces all custom profiles in settings */
    static void migrateOldSettings(QSettings &settings); /* Colors set with old settings keys become a custom profile */
    static QIcon themedIcon(const QString &resource, const QColor &color); /* SVG icon with "currentColor" replaced by color */

private:
    QMap<QString, QColor> colors;
};

#endif // COLORPROFILE_H
