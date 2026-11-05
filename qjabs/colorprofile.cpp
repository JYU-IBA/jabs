#include "colorprofile.h"
#include <QFile>
#include <QGuiApplication>
#include <QPainter>
#include <QPalette>
#include <QPixmap>
#include <QStyleHints>
#include <QSvgRenderer>

const QString ColorProfile::automaticName = QStringLiteral("Automatic");
const QString ColorProfile::lightName = QStringLiteral("Light");
const QString ColorProfile::darkName = QStringLiteral("Dark");

static const QColor reactionColorsDefault[ColorProfile::n_reaction_colors] = { /* Paired colors (ColorBrewer), readable on light and dark backgrounds */
    QColor(0xfb, 0x9a, 0x99),
    QColor(0xe3, 0x1a, 0x1c),
    QColor(0xfd, 0xbf, 0x6f),
    QColor(0xff, 0x7f, 0x00),
    QColor(0xca, 0xb2, 0xd6),
    QColor(0x6a, 0x3d, 0x9a),
    QColor(0xa6, 0xce, 0xe3),
    QColor(0x1f, 0x78, 0xb4),
    QColor(0xb2, 0xdf, 0x8a),
    QColor(0x33, 0xa0, 0x2c),
};

ColorProfile::ColorProfile(const QString &name, bool builtin) : name(name), builtin(builtin)
{
}

const QVector<ColorRole> &ColorProfile::roles()
{
    static QVector<ColorRole> r;
    if(r.isEmpty()) {
        const QColor sys; /* Invalid, system default */
        r = {
            {"editorBackground", "Background", "Editor", QColor("#ffffff"), QColor("#1e1f22"), true},
            {"editorText", "Text", "Editor", QColor("#000000"), QColor("#d8d8d8"), true},
            {"editorComment", "Comment", "Editor", QColor("#6a7a8a"), QColor("#8c9bab"), false},
            {"editorCommand", "Command", "Editor", QColor("#2a62c9"), QColor("#6495ed"), false},
            {"editorVariable", "Variable", "Editor", QColor("#b03a48"), QColor("#f08080"), false},
            {"messageBackground", "Background", "Messages", QColor("#ffffff"), QColor("#1e1f22"), true},
            {"messageText", "Text", "Messages", QColor("#000000"), QColor("#d8d8d8"), true},
            {"messageError", "Error", "Messages", QColor("#cc0000"), QColor("#ff6b6b"), false},
            {"messageWarning", "Warning", "Messages", QColor("#9a6a00"), QColor("#e5c07b"), false},
            {"messageDebug", "Debug", "Messages", QColor("#808080"), QColor("#9a9a9a"), false},
            {"plotBackground", "Background", "Plot", QColor("#ffffff"), QColor("#1e1f22"), false},
            {"plotForeground", "Axes and text", "Plot", QColor("#000000"), QColor("#d8d8d8"), false},
            {"plotGrid", "Grid", "Plot", QColor("#c8c8c8"), QColor("#4a4a4a"), false},
            {"plotLegendBackground", "Legend background", "Plot", QColor(255, 255, 255, 127), QColor(30, 31, 34, 170), false},
            {"plotExperimental", "Experimental", "Plot", QColor("#000000"), QColor("#f0f0f0"), false},
            {"plotSimulated", "Simulated", "Plot", QColor("#0000ff"), QColor("#5b9bff"), false},
            {"plotReference", "Reference", "Plot", QColor("#808080"), QColor("#9a9a9a"), false},
            {"plotConfidence", "Confidence limit", "Plot", QColor(0, 0, 0, 50), QColor(255, 255, 255, 60), false},
        };
        for(int i = 0; i < n_reaction_colors; i++) {
            r.append({QString("plotReaction%1").arg(i + 1), QString("Reaction %1").arg(i + 1), "Plot", reactionColorsDefault[i], reactionColorsDefault[i], false});
        }
        r.append({"iconColor", "Icons", "Toolbar", sys, sys, true});
    }
    return r;
}

QColor ColorProfile::color(const QString &key) const
{
    return colors.value(key, QColor());
}

QColor ColorProfile::color(const QString &key, const QColor &systemColor) const
{
    QColor c = color(key);
    return c.isValid() ? c : systemColor;
}

void ColorProfile::setColor(const QString &key, const QColor &color)
{
    colors[key] = color;
}

Qt::ColorScheme ColorProfile::colorScheme() const
{
    if(name == automaticName) {
        return Qt::ColorScheme::Unknown;
    }
    if(name == lightName) {
        return Qt::ColorScheme::Light;
    }
    if(name == darkName) {
        return Qt::ColorScheme::Dark;
    }
    const QColor bg = color("editorBackground");
    if(!bg.isValid() || !color("messageBackground").isValid()) { /* Custom profile using system backgrounds */
        return Qt::ColorScheme::Unknown;
    }
    return bg.lightness() < 128 ? Qt::ColorScheme::Dark : Qt::ColorScheme::Light;
}

QColor ColorProfile::reactionColor(int index) const
{
    return color(QString("plotReaction%1").arg(index % n_reaction_colors + 1), reactionColorsDefault[index % n_reaction_colors]);
}

ColorProfile ColorProfile::light()
{
    ColorProfile p(lightName, true);
    for(const ColorRole &role : roles()) {
        p.colors[role.key] = role.light;
    }
    return p;
}

ColorProfile ColorProfile::dark()
{
    ColorProfile p(darkName, true);
    for(const ColorRole &role : roles()) {
        p.colors[role.key] = role.dark;
    }
    return p;
}

ColorProfile ColorProfile::automatic()
{ /* Light or dark profile, following the OS. Roles that can be system colors are system colors (e.g. backgrounds), so that they follow the OS. */
    ColorProfile p = systemIsDark() ? dark() : light();
    p.name = automaticName;
    for(const ColorRole &role : roles()) {
        if(role.systemAllowed) {
            p.colors[role.key] = QColor();
        }
    }
    return p;
}

bool ColorProfile::systemIsDark()
{ /* Decided from the palette actually in use, so that profile colors match "System" colors (from the same palette).
   * The color scheme hint (QStyleHints) can be out of sync with the palette while the OS setting is changing. */
    return QGuiApplication::palette().color(QPalette::Window).lightness() < 128;
}

bool ColorProfile::isBuiltinName(const QString &name)
{
    return name == automaticName || name == lightName || name == darkName;
}

QStringList ColorProfile::customProfileNames(QSettings &settings)
{
    settings.beginGroup("colorProfiles");
    QStringList names = settings.childGroups();
    settings.endGroup();
    names.removeIf([](const QString &n) {return isBuiltinName(n);});
    return names;
}

ColorProfile ColorProfile::load(QSettings &settings, const QString &name)
{
    if(name == lightName) {
        return light();
    }
    if(name == darkName) {
        return dark();
    }
    if(name.isEmpty() || name == automaticName || !customProfileNames(settings).contains(name)) {
        return automatic();
    }
    ColorProfile p = light(); /* Roles missing from settings get light defaults */
    p.name = name;
    p.builtin = false;
    settings.beginGroup("colorProfiles");
    settings.beginGroup(name);
    for(const ColorRole &role : roles()) {
        if(settings.contains(role.key)) {
            QString s = settings.value(role.key).toString();
            p.colors[role.key] = (s.isEmpty() && role.systemAllowed) ? QColor() : QColor::fromString(s);
        }
    }
    settings.endGroup();
    settings.endGroup();
    return p;
}

ColorProfile ColorProfile::active(QSettings &settings)
{
    return load(settings, settings.value("colorProfile", automaticName).toString());
}

void ColorProfile::saveCustomProfiles(QSettings &settings, const QList<ColorProfile> &profiles)
{
    settings.remove("colorProfiles");
    settings.beginGroup("colorProfiles");
    for(const ColorProfile &p : profiles) {
        if(p.builtin || isBuiltinName(p.name)) {
            continue;
        }
        settings.beginGroup(p.name);
        for(const ColorRole &role : roles()) {
            QColor c = p.color(role.key);
            settings.setValue(role.key, c.isValid() ? c.name(c.alpha() == 255 ? QColor::HexRgb : QColor::HexArgb) : QString());
        }
        settings.endGroup();
    }
    settings.endGroup();
}

void ColorProfile::migrateOldSettings(QSettings &settings)
{
    if(settings.value("colorProfilesMigrated", false).toBool()) {
        return;
    }
    settings.setValue("colorProfilesMigrated", true);
    const QStringList oldKeys = {"messageColor", "commentColor", "commandColor", "variableColor"};
    const QStringList newKeys = {"messageText", "editorComment", "editorCommand", "editorVariable"};
    bool found = false;
    ColorProfile p = systemIsDark() ? dark() : light();
    p.name = "Custom";
    p.builtin = false;
    for(int i = 0; i < oldKeys.size(); i++) {
        QString s = settings.value(oldKeys[i]).toString();
        if(!s.isEmpty() && QColor::fromString(s).isValid()) {
            p.colors[newKeys[i]] = QColor::fromString(s);
            found = true;
        }
    }
    if(!found || customProfileNames(settings).contains(p.name)) {
        return;
    }
    settings.beginGroup("colorProfiles");
    settings.beginGroup(p.name);
    for(const ColorRole &role : roles()) {
        QColor c = p.color(role.key);
        settings.setValue(role.key, c.isValid() ? c.name(c.alpha() == 255 ? QColor::HexRgb : QColor::HexArgb) : QString());
    }
    settings.endGroup();
    settings.endGroup();
    settings.setValue("colorProfile", p.name);
}

QIcon ColorProfile::themedIcon(const QString &resource, const QColor &color)
{
    QFile file(resource);
    if(!file.open(QIODevice::ReadOnly)) {
        return QIcon(resource);
    }
    QByteArray svg = file.readAll();
    svg.replace("currentColor", color.name(QColor::HexRgb).toUtf8());
    QSvgRenderer renderer(svg);
    if(!renderer.isValid()) {
        return QIcon(resource);
    }
    QIcon icon;
    for(int size : {16, 24, 32, 48, 64, 128}) {
        QPixmap pixmap(size, size);
        pixmap.fill(Qt::transparent);
        QPainter painter(&pixmap);
        renderer.render(&painter);
        painter.end();
        icon.addPixmap(pixmap);
    }
    return icon;
}
