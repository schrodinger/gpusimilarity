"""
This is a sample HTTP server interface with the FastSim backend,
which takes fingerprints as a JSON and returns results in JSON form.
"""

import os
import subprocess
import time

from PyQt5 import QtCore, QtNetwork
from PyQt5 import hunt

from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn

import cgi
import json

import fastsim_utils
from fastsim_utils import RETURN_COUNT

SCRIPT_DIR = os.path.split(__file__)[0]
BITCOUNT = 1024
sockets = {}


class ThreadedHTTPServer(ThreadingMixIn, HTTPServer):
    """Use threads to handle requests"""


class FastSimHandler(BaseHTTPRequestHandler):

    """
    Retrieve the smiles passed into the form and the reults from the backend
    """
    def get_data(self, socket_name):
        form = cgi.FieldStorage(
            fp=self.rfile,
            headers=self.headers,
            environ={'REQUEST_METHOD': 'POST',
                     'CONTENT_TYPE': self.headers['Content-Type'],
                     })

        search_smiles = form["smiles"].value.strip()
        fp_binary = fastsim_utils.smiles_to_fingerprint_bin(search_smiles)
        fp_qba = QtCore.QByteArray(fp_binary)

        socket = sockets[socket_name]
        socket.write(fp_qba)
        socket.flush()
        socket.waitForReadyRead(5000)

        output_qba = socket.readAll()
        return search_smiles, output_qba

    def do_GET(self):
        self.send_error(404, 'Server unavailable.')

    def do_POST(self):
        if not self.path.startswith("/similarity_search_json"):
            return
        dbname = self.path.replace("/similarity_search_json_", "")
        self.send_error(404, 'DB %s not found in options: %s' % (dbname, sockets.keys())) #noqa
        return

        _, output_qba = self.get_data(dbname)
        self.do_json_POST(output_qba)

    def results2json(self, output_qba):
        data_reader = QtCore.QDataStream(output_qba)
        smiles, ids, scores = [], [], []
        for i in range(RETURN_COUNT):
            smiles.append(data_reader.readString().decode("utf-8"))
        for i in range(RETURN_COUNT):
            ids.append(data_reader.readString().decode("utf-8"))
        for i in range(RETURN_COUNT):
            scores.append(data_reader.readFloat())
        results = []
        for cid, smi, score in zip(ids, smiles, scores):
            results.append([cid, smi, score])
        return results

    def do_json_POST(self, output_qba):
        self.send_response(200)
        self.send_header('Content-type', 'text/json')
        self.end_headers()

        results_json = self.results2json(output_qba)
        self.wfile.write(json.dumps(results_json).encode('utf8'))


class FastSimHTTPHandler(FastSimHandler):
    def _generateImageIfNeeded(self, fullpath):
        if (self.path.startswith('smiles_') and not os.path.exists(fullpath)):
            safe_smi = self.path.replace('smiles_', '').replace('.png', '')
            smi = safe_smi.replace('_-1-_', '/').replace(
                '_-2-_', '\\').replace('_-3-_', '#')
            fastsim_utils.smiles_to_image_file(smi, fullpath)

    def write_results_html(self, search_smiles, smiles, scores, ids):
        for smi, score, cid in zip(smiles, scores, ids):
            id_html = cid
            if cid.startswith('ZINC'):
                id_html = \
                    "<a href=http://zinc.docking.org/substance/{0}>{1}</a>".format(cid[4:], cid) #noqa

            safe_smi = smi.replace('/', '_-1-_').replace(
                '\\', '_-2-_').replace('#', '_-3-_')
            self.wfile.write("""<img src='smiles_{0}.png'>
                                <img src='smiles_{1}.png'>
                                <table><tr><td>{3}: {1}</td></tr>
                                <tr><td>{2}</td></tr><td></table>""".format(
                                    search_smiles, safe_smi, score,
                                    id_html).encode())

    def do_GET(self):
        if self.path == "/":
            self.path = "/index.html"
        if self.path.startswith("/"):
            self.path = self.path[1:]

        try:
            sendReply = False
            fullpath = os.path.join(SCRIPT_DIR, self.path)
            if self.path.endswith(".html"):
                mimetype = 'text/html'
                sendReply = True
            if self.path.endswith(".png"):
                mimetype = 'image/png'
                sendReply = True
                self._generateImageIfNeeded(fullpath)

            if sendReply is True:
                f = open(fullpath, 'rb+')
                self.send_response(200)
                self.send_header('Content-type', mimetype)
                self.end_headers()
                self.wfile.write(f.read())
                f.close()
            return

        except IOError:
            self.send_error(404, 'File Not Found: %s' % self.path)

    def do_POST(self):
        if not self.path.startswith("/similarity_search"):
            return

        dbname = self.path.replace("/similarity_search_", "")
        dbname = dbname.replace("json_", "")  # Strip json if in endpoint
        if dbname not in sockets:
            self.send_error(404, 'DB %s not found in options: %s' % (dbname, sockets.keys())) #noqa
            return

        search_smiles, output_qba = self.get_data(dbname)
        if self.path.startswith("/similarity_search_json"):
            self.do_json_POST(output_qba)
        elif self.path.startswith("/similarity_search"):
            self.do_html_POST(output_qba, search_smiles)

    def deserialize_results(self, output_qba):
            data_reader = QtCore.QDataStream(output_qba)
            smiles, ids, scores = [], [], []
            for i in range(RETURN_COUNT):
                smiles.append(data_reader.readString().decode("utf-8"))
            for i in range(RETURN_COUNT):
                ids.append(data_reader.readString().decode("utf-8"))
            for i in range(RETURN_COUNT):
                scores.append(data_reader.readFloat())
            return smiles, ids, scores

    def do_html_POST(self, output_qba, search_smiles):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

        smiles, ids, scores = self.deserialize_results(output_qba)
        f = open(os.path.join(SCRIPT_DIR, "index.html"), 'rb+')
        self.wfile.write(f.read())
        self.write_results_html(search_smiles, smiles, scores, ids)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Sample FastSim Server - "
            "run an HTTP server that loads fingerprint data onto GPU and " #noqa
            "responds to queries to find most similar fingperints.") #noqa
    parser.add_argument('dbnames', help=".fsim files containing fingerprint "
                        "data to be searched", nargs='*')
    parser.add_argument('--hostname', default="localhost",
                        help="Hostname to run on")
    parser.add_argument('--port', default=8080, type=int,
                        help="Port to run on")
    parser.add_argument('--http_interface', action='store_true',
                        help="Start HTTP server for debugging, not secure enough for production machine") #noqa
    return parser.parse_args()


def main():

    global sockets
    args = parse_args()

    # Try to connect to the GPU backend
    app = QtCore.QCoreApplication([])

    procs = []
    for dbname in args.dbnames:
        # Start the GPU backend
        combi_bin_dir = hunt('combiglide', 'exec')
        fastsim_exec = os.path.join(combi_bin_dir, 'fastsim_server')
        procs.append(subprocess.Popen([fastsim_exec, dbname]))
        socket = QtNetwork.QLocalSocket(app)
        dbname_noext = os.path.splitext(dbname)[0]
        sockets[dbname_noext] = socket
        while not socket.isValid():
            socket_name = os.path.splitext(os.path.basename(dbname))[0]
            socket.connectToServer(socket_name)
            time.sleep(0.3)

    if args.http_interface:
        handler = FastSimHTTPHandler
    else:
        handler = FastSimHandler
    server = ThreadedHTTPServer((args.hostname, args.port), handler)
    print("Running HTTP server...")
    server.serve_forever()
    for proc in procs:
        proc.kill()


if __name__ == '__main__':
    main()
