from PyQt5 import QtCore, QtNetwork

import random
from gpusim_utils import smiles_to_fingerprint_bin


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Sample GPUSim Server - "
            "run an HTTP server that loads fingerprint data onto GPU and " #noqa
            "responds to queries to find most similar fingperints.") #noqa
    parser.add_argument('dbname', help=".fsim file containing fingerprint "
                        "data to be searched")
    parser.add_argument('dbkey', default="", help="Key for fsim file")
    return parser.parse_args()


def main():
    args = parse_args()
    app = QtCore.QCoreApplication([])

    socket = QtNetwork.QLocalSocket(app)
    smiles = input("Smiles: ")
    dbcount = 1
    dbname = args.dbname
    dbkey = args.dbkey
    socket.connectToServer('gpusimilarity')

    while smiles and smiles.lower() not in ('quit', 'exit'):
        return_count = 20
        similarity_cutoff = 0

        fp_binary, _ = smiles_to_fingerprint_bin(smiles)
        fp_qba = QtCore.QByteArray(fp_binary)

        output_qba = QtCore.QByteArray()
        output_qds = QtCore.QDataStream(output_qba, QtCore.QIODevice.WriteOnly)

        output_qds.writeInt(dbcount)
        output_qds.writeString(dbname.encode())
        output_qds.writeString(dbkey.encode())

        request_num = random.randint(0, 2**31)
        output_qds.writeInt(request_num)
        output_qds.writeInt(return_count)
        output_qds.writeFloat(similarity_cutoff)
        output_qds << fp_qba

        socket.write(output_qba)
        socket.flush()
        socket.waitForReadyRead(30000)
        output_qba = socket.readAll()

        smiles = []
        scores = []
        ids = []

        data_reader = QtCore.QDataStream(output_qba)
        returned_request = data_reader.readInt()
        if request_num != returned_request:
            raise RuntimeError("Incorrect result ID returned!")

        return_count = data_reader.readInt()
        approximate_matches = data_reader.readUInt64()

        for i in range(return_count):
            smiles.append(data_reader.readString())
        for i in range(return_count):
            ids.append(data_reader.readString())
        for i in range(return_count):
            scores.append(data_reader.readFloat())

        print("Approximate total matches: {0}, returning {1}".format(
            approximate_matches, return_count))
        for cid, smi, score in zip(ids, smiles, scores):
            print("{0} {1}: {2}".format(cid, smi, score))
        smiles = input("Smiles: ")


if __name__ == '__main__':
    main()
