#include <QCoreApplication>
#include <QFileInfo>
#include <iostream>

#include "gpusim.h"

using gpusim::GPUSimServer;
using gpusim::Fingerprint;

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    bool file_exists = QFileInfo(argv[1]).exists();
    if(!file_exists) {
        std::cout << "File: \"" << argv[1] << "\" not found." << std::endl;
    }

    if(app.arguments().contains("--help") || app.arguments().contains("-h") ||
            !file_exists) {
        std::cout << "Arg parsing is only done in a reasonable way in"
            " the python gpusim_server.py.  Handling here is very error prone"
            " and not intended for direct use." << std::endl;
        return 1;
    }
    GPUSimServer gpusim(argv[1]);
    if(app.arguments().contains("--cpu_only")) {
        gpusim.setUseGPU(false);
    }
    app.exec();

    return 0;
};
