"""
This is a sample HTTP server interface with the GPUSim backend,
which takes fingerprints as a JSON and returns results in JSON form.

THIS RUNS ON PORT 80 AND IS NOT SECURE.  For production use you should wrap
using nginx or equivalent and only use https externally.
"""

import os
import random
import subprocess
import sys
import time
import tempfile

from PyQt5 import QtCore, QtNetwork

from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn

import cgi
import json

import gpusim_utils

SCRIPT_DIR = os.path.split(__file__)[0]
BITCOUNT = 1024

socket = None

# Make sure there's only ever a single search at a time
search_mutex = QtCore.QMutex()

try:
    from gpusim_server_loc import GPUSIM_EXEC  # Used in schrodinger env
except ImportError:
    GPUSIM_EXEC = './gpusimserver'


class ThreadedHTTPServer(ThreadingMixIn, HTTPServer):
    """Use threads to handle requests"""


class GPUSimHandler(BaseHTTPRequestHandler):

    """
    Retrieve the smiles passed into the form and the results from the backend
    """


    def get_posted_data(self):
        form = cgi.FieldStorage(
            fp=self.rfile,
            headers=self.headers,
            environ={'REQUEST_METHOD': 'POST',
                     'CONTENT_TYPE': self.headers['Content-Type'],
                     })

        dbnames = form["dbnames"].value.split(',')
        try:
            dbkeys = form["dbkeys"].value.split(',')
        except KeyError:  # Handle empty string not passed by HTML post
            dbkeys = [""]
        if len(dbnames) != len(dbkeys):
            raise RuntimeError("Need key for each database.")

        return (form["smiles"].value.strip(), int(form["return_count"].value),
                float(form["similarity_cutoff"].value), dbnames, dbkeys)

    def get_data(self, dbnames, dbkeys, src_smiles, return_count,
                 similarity_cutoff, request_num):
        global socket
        fp_binary, canon_smile = gpusim_utils.smiles_to_fingerprint_bin(src_smiles)
        fp_qba = QtCore.QByteArray(fp_binary)

        output_qba = QtCore.QByteArray()
        output_qds = QtCore.QDataStream(output_qba, QtCore.QIODevice.WriteOnly)

        output_qds.writeInt(len(dbnames))
        for name, key in zip(dbnames, dbkeys):
            output_qds.writeString(name.encode())
            output_qds.writeString(key.encode())

        output_qds.writeInt(request_num)
        output_qds.writeInt(return_count)
        output_qds.writeFloat(similarity_cutoff)
        output_qds << fp_qba

        socket.write(output_qba)
        socket.flush()
        socket.waitForReadyRead(30000)

        output_qba = socket.readAll()
        return output_qba

    def do_GET(self):
        self.send_error(404, 'Server unavailable.')

    def search_for_results(self):
        global search_mutex
        search_mutex.lock()
        try:
            src_smiles, return_count, similarity_cutoff, dbnames, dbkeys = \
                self.get_posted_data()
            request_num = random.randint(0, 2**31)
            print("Processing request {0}".format(request_num),
                  file=sys.stderr) #noqa
            output_qba = self.get_data(dbnames, dbkeys, src_smiles,
                                       return_count, similarity_cutoff,
                                       request_num)
            try:
                return_count, smiles, ids, scores = self.deserialize_results(
                    request_num, output_qba)
            except RuntimeError:
                self.flush_socket()
                raise

            return smiles, ids, scores, src_smiles
        finally:
            search_mutex.unlock()

    def flush_socket(self):
        global socket
        while not socket.atEnd():
            socket.readAll()

    def do_POST(self):
        if not self.path.startswith("/similarity_search_json"):
            return
        try:
            smiles, ids, scores, _ = self.search_for_results()
        except ValueError:
            return
        self.do_json_POST(smiles, ids, scores)

    def deserialize_results(self, request_num, output_qba):
        data_reader = QtCore.QDataStream(output_qba)
        returned_request = data_reader.readInt()
        if request_num != returned_request:
            raise RuntimeError("Incorrect result ID returned!")
        return_count = data_reader.readInt()
        smiles, ids, scores = [], [], []
        for i in range(return_count):
            smiles.append(data_reader.readString().decode("utf-8"))
        for i in range(return_count):
            ids.append(data_reader.readString().decode("utf-8"))
        for i in range(return_count):
            scores.append(data_reader.readFloat())
        return return_count, smiles, ids, scores

    def results2json(self, smiles, ids, scores):
        results = []
        for cid, smi, score in zip(ids, smiles, scores):
            results.append([cid, smi, score])
        return results

    def do_json_POST(self, smiles, ids, scores):
        self.send_response(200)
        self.send_header('Content-type', 'text/json')
        self.end_headers()

        results_json = self.results2json(smiles, ids, scores)
        self.wfile.write(json.dumps(results_json).encode('utf8'))


class GPUSimHTTPHandler(GPUSimHandler):
    _tmp_dir = tempfile.TemporaryDirectory()

    def _generateImageIfNeeded(self, fullpath):
        if (self.path.startswith('smiles_') and not os.path.exists(fullpath)):
            safe_smi = self.path.replace('smiles_', '').replace('.png', '')
            smi = safe_smi.replace('_-1-_', '/').replace(
                '_-2-_', '\\').replace('_-3-_', '#')
            gpusim_utils.smiles_to_image_file(smi, fullpath)

    def write_results_html(self, src_smiles, smiles, scores, ids):
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
                                    src_smiles, safe_smi, score,
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
                fullpath = os.path.join(self._tmp_dir.name, self.path)
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

        try:
            smiles, ids, scores, src_smiles = self.search_for_results()
        except ValueError:
            return
        if self.path.startswith("/similarity_search_json"):
            self.do_json_POST(smiles, ids, scores)
        elif self.path.startswith("/similarity_search"):
            self.do_html_POST(smiles, ids, scores, src_smiles)

    def do_html_POST(self, smiles, ids, scores, src_smiles):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

        f = open(os.path.join(SCRIPT_DIR, "index.html"), 'rb+')
        self.wfile.write(f.read())
        self.write_results_html(src_smiles, smiles, scores, ids)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Sample GPUSim Server - "
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
    parser.add_argument('--cpu_only', action='store_true',
                        help="Search the database on the CPU, not the GPU (slow)") #noqa
    parser.add_argument('--gpu_bitcount', default='0',
                        help="Provide the maximum bitcount for fingerprints on GPU") #noqa
    parser.add_argument('--debug', action='store_true', help="Run the backend inside GDB") #noqa
    return parser.parse_args()


def setup_socket(app):
    global socket

    socket = QtNetwork.QLocalSocket(app)
    while not socket.isValid():
        socket_name = 'gpusimilarity'
        socket.connectToServer(socket_name)
        time.sleep(0.3)


def main():

    args = parse_args()

    # Try to connect to the GPU backend
    app = QtCore.QCoreApplication([])

    # Start the GPU backend
    cmdline = [GPUSIM_EXEC]
    if args.cpu_only:
        cmdline.append('--cpu_only')
    cmdline += ['--gpu_bitcount', args.gpu_bitcount]
    cmdline += args.dbnames
    backend_proc = subprocess.Popen(cmdline)
    setup_socket(app)

    if args.http_interface:
        handler = GPUSimHTTPHandler
    else:
        handler = GPUSimHandler
    server = ThreadedHTTPServer((args.hostname, args.port), handler)
    print("Running HTTP server...", file=sys.stderr)
    try:
        server.serve_forever()
    finally:
        backend_proc.kill()


if __name__ == '__main__':
    main()
