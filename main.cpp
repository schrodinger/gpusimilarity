#include <QCommandLineParser>
#include <QCoreApplication>
#include <QDebug>
#include <QFileInfo>
#include <iostream>

#include "gpusim.h"

using gpusim::Fingerprint;
using gpusim::GPUSimServer;

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("GPUSimilarity");
    QCommandLineParser parser;
    parser.setApplicationDescription("GPUSimilarity Backend:  Not meant to be "
                                     "called directly.");
    const QCommandLineOption helpOption = parser.addHelpOption();

    QCommandLineOption cpuOnlyOption("cpu_only",
                                     "Perform searches only on CPU");
    parser.addOption(cpuOnlyOption);

    QCommandLineOption gpuBitcountOption(
        "gpu_bitcount", "Define the fingerprint bitcount on the GPU",
        "Bitcount", "0");
    parser.addOption(gpuBitcountOption);

    if (!parser.parse(QCoreApplication::arguments())) {
        qDebug() << parser.errorText();
        return 1;
    }

    if (parser.isSet(helpOption)) {
        qDebug() << "Arg parsing is only done in a reasonable way in"
                    " the python gpusim_server.py.  Handling here is very "
                    "error prone"
                    " and not intended for direct use.";
        return 1;
    }

    bool ok = false;
    auto gpu_bitcount = parser.value(gpuBitcountOption).toInt(&ok);
    if (!ok) {
        qDebug() << "GPU Bitcount must be an integer";
        return 1;
    }

    if (parser.isSet(cpuOnlyOption) && gpu_bitcount != 0) {
        qDebug() << "--cpu_only and --gpu_bitcount are incompatible options";
        return 1;
    }

    auto db_fnames = parser.positionalArguments();
    for (auto filename : db_fnames) {
        bool file_exists = QFileInfo(filename).exists();
        if (!file_exists) {
            qDebug() << "File: \"" << qPrintable(filename) << "\" not found.";
            return 1;
        }
    }

    GPUSimServer gpusim(db_fnames, gpu_bitcount);
    gpusim.setUseGPU(!parser.isSet(cpuOnlyOption));

    app.exec();

    return 0;
};
