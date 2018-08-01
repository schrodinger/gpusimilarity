#include <QCoreApplication>
#include <QFileInfo>
#include <iostream>

#include "gpusim.h"

using gpusim::GPUSimServer;
using gpusim::Fingerprint;

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    auto args = app.arguments();
    bool use_gpu = true;

    if(args.contains("--cpu_only")) {
        use_gpu = false;
    }

    int first_file_idx = 1;
    while(args[first_file_idx].startsWith("--")) first_file_idx++;

    bool file_exists = QFileInfo(args[first_file_idx]).exists();
    if(!file_exists) {
        std::cout << "File: \"" << qPrintable(args[first_file_idx]) << "\" not found." << std::endl;
    }

    if(args.contains("--help") || args.contains("-h") ||
            !file_exists) {
        std::cout << "Arg parsing is only done in a reasonable way in"
            " the python gpusim_server.py.  Handling here is very error prone"
            " and not intended for direct use." << std::endl;
        return 1;
    }
    GPUSimServer gpusim(args.mid(first_file_idx));
    gpusim.setUseGPU(use_gpu);

    app.exec();

    return 0;
};
