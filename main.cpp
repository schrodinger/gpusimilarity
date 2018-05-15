#include <QCoreApplication>

#include "fastsim.h"

using fastsim::FastSimServer;
using fastsim::Fingerprint;

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    FastSimServer fastsim(argv[1]);
    if(app.arguments().contains("--cpu_only")) {
        fastsim.setUseGPU(false);
    }
    app.exec();

    return 0;
};
