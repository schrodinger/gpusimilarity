#pragma once

#if QT_VERSION < QT_VERSION_CHECK(5, 5, 0)
#define qInfo qDebug
#define qWarning qDebug
#define qFatal qDebug
#define qCritical qDebug
#endif
