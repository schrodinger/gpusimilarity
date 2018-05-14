"""
fastsim_createdb.py
Use this to create a serialized data file to be used with the fastsim backend.

NOTE:  GPGPU backend requires fingerprint size to be sizeof(int) divisible
"""


from PyQt5 import QtCore
import gzip

import fastsim_utils


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='A toy!')
    parser.add_argument('inputfile')
    parser.add_argument('outputfile')
    return parser.parse_args()


def main():
    args = parse_args()
    qf = QtCore.QFile(args.outputfile)
    qf.open(QtCore.QIODevice.WriteOnly)

    count = 0

    print("Processing file {0}".format(args.inputfile))
    input_fhandle = gzip.open(args.inputfile, 'rb')
    print("Processing Smiles...")

    smi_byte_data = QtCore.QByteArray()
    smi_qds = QtCore.QDataStream(smi_byte_data, QtCore.QIODevice.WriteOnly)

    id_byte_data = QtCore.QByteArray()
    id_qds = QtCore.QDataStream(id_byte_data, QtCore.QIODevice.WriteOnly)

    fp_byte_data = QtCore.QByteArray()
    fp_qds = QtCore.QDataStream(fp_byte_data, QtCore.QIODevice.WriteOnly)

    print("Reading lines...")
    read_bytes = 10000000
    count = 0
    lines = input_fhandle.readlines(read_bytes)
    print(len(lines))
    while lines != []:

        rows = fastsim_utils.split_lines_add_fp(lines)
        filtered_rows = [row for row in rows if row is not None]
        count += len(filtered_rows)
        for row in filtered_rows:
            if row is None:
                continue
            smi_qds.writeString(row[0])
            id_qds.writeString(row[1])
            fp_qds.writeRawData(row[2])

        print("Processed {0} rows".format(count))
        lines = input_fhandle.readlines(read_bytes)

    qds = QtCore.QDataStream(qf)
    # Set version so that files will be usable cross-release
    qds.setVersion(QtCore.QDataStream.Qt_5_2)

    size = QtCore.QSize(fastsim_utils.BITCOUNT, count)
    qds << size
    qds << fp_byte_data
    qds << smi_byte_data
    qds << id_byte_data

    qf.close()


if __name__ == "__main__":
    main()
