"""
gpusim_createdb.py
Use this to create a serialized data file to be used with the gpusim backend.

NOTE:  GPGPU backend requires fingerprint size to be sizeof(int) divisible
"""

from PyQt5 import QtCore
import gzip

import gpusim_utils

DATA_VERSION = 2


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Create a GPUSimilarity '
                                     'Binary FingerprintDB')
    parser.add_argument('inputfile')
    parser.add_argument('outputfile')
    parser.add_argument('--trustSmiles', action='store_true', default=False)
    return parser.parse_args()


try:
    import ipyparallel as ipp
    rc = ipp.Client()
    dview = rc[:]
except ImportError:
    dview = None


class FPData:
    """
    A class to store all the serialized Fingerprint data, as well as smiles and
    ID data, of the database
    """

    def __init__(self):
        self.smi_byte_data = [QtCore.QByteArray()]
        self.smi_qds = QtCore.QDataStream(self.smi_byte_data[0],
                                          QtCore.QIODevice.WriteOnly)

        self.id_byte_data = [QtCore.QByteArray()]
        self.id_qds = QtCore.QDataStream(self.id_byte_data[0],
                                         QtCore.QIODevice.WriteOnly)

        self.fp_byte_data = [QtCore.QByteArray()]
        self.fp_qds = QtCore.QDataStream(self.fp_byte_data[0],
                                         QtCore.QIODevice.WriteOnly)

    def checkQBASize(self, qba_list, qds):
        """
        If the current QByteArray has hit its limit, create a new one, and return
        the new QDataStream to that array
        @param qba_list: List of QByteArrays holding relevant data
        @type  qba_list: [QByteArray, ...]
        @param qds: Stream to currently used QByteArray
        @type  qds: QDataStream
        """
        if qba_list[-1].size() < 2**30:
            return qds

        qba_list.append(QtCore.QByteArray())
        return QtCore.QDataStream(qba_list[-1], QtCore.QIODevice.WriteOnly)

    def checkQBASizes(self):
        """ Make sure the QByteArray lists are in correct size/state """
        self.smi_qds = self.checkQBASize(self.smi_byte_data, self.smi_qds)
        self.id_qds = self.checkQBASize(self.id_byte_data, self.id_qds)
        self.fp_qds = self.checkQBASize(self.fp_byte_data, self.fp_qds)

    def storeData(self, row):
        """ A single rows data to be serialized/stored """
        if row is None:
            return
        self.checkQBASizes()
        self.smi_qds.writeString(row[0])
        self.id_qds.writeString(row[1])
        self.fp_qds.writeRawData(row[2])

    def writeData(self, qds):
        """
        Write all the saved serialized data to a QDataStream
        @param qds: Stream output to file
        @type  qds: QDataStream
        """
        for qba_list in [self.fp_byte_data, self.smi_byte_data,
                         self.id_byte_data]:
            qds.writeInt(len(qba_list))
            for byte_data in qba_list:
                qds << byte_data


def main():
    args = parse_args()
    qf = QtCore.QFile(args.outputfile)
    qf.open(QtCore.QIODevice.WriteOnly)

    count = 0

    print("Processing file {0}".format(args.inputfile))
    input_fhandle = gzip.open(args.inputfile, 'rb')
    print("Processing Smiles...")

    fpdata = FPData()

    print("Reading lines...")
    read_bytes = 10000000
    count = 0
    lines = input_fhandle.readlines(read_bytes)
    print(len(lines))
    while lines != []:
        rows = gpusim_utils.split_lines_add_fp(
            lines, dview=dview, trust_smiles=args.trustSmiles)
        filtered_rows = [row for row in rows if row is not None]
        count += len(filtered_rows)
        for row in filtered_rows:
            fpdata.storeData(row)

        print("Processed {0} rows".format(count))
        lines = input_fhandle.readlines(read_bytes)

    qds = QtCore.QDataStream(qf)
    # Set version so that files will be usable cross-release
    qds.setVersion(QtCore.QDataStream.Qt_5_2)

    qds.writeInt(DATA_VERSION)
    qds.writeInt(gpusim_utils.BITCOUNT)
    qds.writeInt(count)
    fpdata.writeData(qds)

    qf.close()


if __name__ == "__main__":
    main()
