import sys
from PyQt5 import QtCore, QtNetwork

from fastsim_utils import smiles_to_fingerprint_bin, RETURN_COUNT


def main():
    app = QtCore.QCoreApplication([])

    socket = QtNetwork.QLocalSocket(app)
    smiles = input("Smiles: ")
    db_name = sys.argv[1]
    socket.connectToServer(db_name)

    while smiles and smiles.lower() not in ('quit', 'exit'):
        fp_binary = smiles_to_fingerprint_bin(smiles)
        fp_qba = QtCore.QByteArray(fp_binary)

        socket.write(fp_qba)
        socket.flush()
        socket.waitForReadyRead(5000)
        output_qba = socket.readAll()

        smiles = []
        scores = []
        ids = []

        data_reader = QtCore.QDataStream(output_qba)

        for i in range(RETURN_COUNT):
            smiles.append(data_reader.readString())
        for i in range(RETURN_COUNT):
            ids.append(data_reader.readString().decode("utf-8"))
        for i in range(RETURN_COUNT):
            scores.append(data_reader.readFloat())

        for cid, smi, score in zip(ids, smiles, scores):
            print("{0} {1}: {2}".format(cid, smi, score))
        smiles = input("Smiles: ")


if __name__ == '__main__':
    main()
