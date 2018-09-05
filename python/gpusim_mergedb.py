"""
gpusim_mergedb.py
Use this to merge fsim files into one, allowing for highly parallel creation

Example command:
python3 gpusmi_mergedb.py in1.fsim in2.fsim in3.fsim -o out.fsim
"""

from PyQt5 import QtCore

from gpusim_createdb import DATABASE_VERSION


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

    smi_byte_data = []
    id_byte_data = []
    fp_byte_data = []

    count = 0
    bitcount = None
    for dbname in args.dbnames:
        dbf = QtCore.QFile(dbname)
        dbf.open(QtCore.QIODevice.ReadOnly)
        dbds = QtCore.QDataStream(dbf)
        version = dbds.readInt()
        if version != DATABASE_VERSION:
            raise RuntimeError("Input database not compatible with this "
                               "version of GPUSim")
        in_bitcount = dbds.readInt()
        in_count = dbds.readInt()
        if bitcount is None:
            bitcount = in_bitcount
        elif bitcount != in_bitcount:
            raise ValueError("Can't mix databases with different "
                             "fingerprint bitcounts")
        for qba_list in [fp_byte_data, smi_byte_data, id_byte_data]:
            list_len = dbds.readInt()
            for i in range(list_len):
                byte_data = QtCore.QByteArray()
                dbds >> byte_data
                qba_list.append(byte_data)

        count += in_count

    print("Writing new database with {0} entries.".format(count))
    qds.writeInt(DATABASE_VERSION)
    qds.writeInt(bitcount)
    qds.writeInt(count)
    for qba_list in [fp_byte_data, smi_byte_data, id_byte_data]:
        qds.writeInt(len(qba_list))
        for byte_data in qba_list:
            qds << byte_data

    qf.close()


if __name__ == "__main__":
    main()
