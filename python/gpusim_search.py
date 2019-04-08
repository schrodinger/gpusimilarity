import sys
from PyQt5 import QtCore, QtNetwork

from gpusim_utils import smiles_to_fingerprint_bin


def main():
    app = QtCore.QCoreApplication([])

    socket = QtNetwork.QLocalSocket(app)
    smiles = input("Smiles: ")
    dbcount = 1
    dbname = sys.argv[1]
    dbkey = sys.argv[2]
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
        return_count = data_reader.readInt()

        for i in range(return_count):
            smiles.append(data_reader.readString())
        for i in range(return_count):
            ids.append(data_reader.readString())
        for i in range(return_count):
            scores.append(data_reader.readFloat())

        for cid, smi, score in zip(ids, smiles, scores):
            print("{0} {1}: {2}".format(cid, smi, score))
        smiles = input("Smiles: ")


if __name__ == '__main__':
    main()
