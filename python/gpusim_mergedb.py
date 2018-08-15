"""
gpusim_mergedb.py
Use this to merge fsim files into one, allowing for highly parallel creation

Example command:
python3 gpusmi_mergedb.py in1.fsim in2.fsim in3.fsim -o out.fsim
"""

from PyQt5 import QtCore

import gpusim_utils


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Merge GPUSimilarity'
                                     ' Binary FingerprintDBs')
    parser.add_argument('--outputfile', '-o', help='Filename to merge fsim '
                        'files into')
    parser.add_argument('dbnames', help=".fsim files containing fingerprint "
                        "data to be merged", nargs='*')
    return parser.parse_args()


def main():
    args = parse_args()
    qf = QtCore.QFile(args.outputfile)
    qf.open(QtCore.QIODevice.WriteOnly)

    qds = QtCore.QDataStream(qf)
    # Set version so that files will be usable cross-release
    qds.setVersion(QtCore.QDataStream.Qt_5_2)

    smi_byte_data = QtCore.QByteArray()
    id_byte_data = QtCore.QByteArray()
    fp_byte_data = QtCore.QByteArray()

    count = 0
    bitcount = None
    for dbname in args.dbnames:
        dbf = QtCore.QFile(dbname)
        dbf.open(QtCore.QIODevice.ReadOnly)
        dbds = QtCore.QDataStream(dbf)
        size = QtCore.QSize(0, 0)
        dbds >> size
        if bitcount is None:
            bitcount = size.width()
        elif bitcount != size.width():
            raise ValueError("Can't mix databases with different "
                             "fingerprint bitcounts")
        local_smi_byte_data = QtCore.QByteArray()
        local_id_byte_data = QtCore.QByteArray()
        local_fp_byte_data = QtCore.QByteArray()
        dbds >> local_fp_byte_data
        dbds >> local_smi_byte_data
        dbds >> local_id_byte_data
        fp_byte_data.append(local_fp_byte_data)
        smi_byte_data.append(local_smi_byte_data)
        id_byte_data.append(local_id_byte_data)
        count += size.height()

    print("Writing new database with {0} entries.".format(count))
    size = QtCore.QSize(gpusim_utils.BITCOUNT, count)
    qds << size
    qds << fp_byte_data
    qds << smi_byte_data
    qds << id_byte_data

    qf.close()


if __name__ == "__main__":
    main()
