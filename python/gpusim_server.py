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
import tempfile
import time

from PyQt5 import QtCore, QtNetwork

from collections import namedtuple
from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn

import cgi
import json

import gpusim_utils

BITCOUNT = 1024
MAX_RETRY = 3  # Used for sporadic QSharedMemory::initKey errors
SCRIPT_DIR = os.path.dirname(__file__)
SERVER_RESULT_TIMEOUT = 5  # time.time() is in seconds
SERVER_ERROR_RESULT = (1, ['CC'], ['SERVER_ERROR_ON_SEARCH'], [0])

socket = None

# Make sure there's only ever a single search at a time
search_mutex = QtCore.QMutex()

QueryParams = namedtuple('QueryParams', 'dbnames dbkeys version '
                         'smiles return_count similarity_cutoff')

try:
    from gpusim_server_loc import GPUSIM_EXEC  # Used in schrodinger env
except ImportError:
    script_path = os.path.dirname(__file__)
    GPUSIM_EXEC = os.path.join(script_path, '..', 'gpusimserver')


class ThreadedHTTPServer(ThreadingMixIn, HTTPServer):
    """Use threads to handle requests"""


class GPUSimHandler(BaseHTTPRequestHandler):

    """
    Retrieve the smiles passed into the form and the results from the backend
    """

    def get_posted_query(self):
        form = cgi.FieldStorage(
            fp=self.rfile,
            headers=self.headers,
            environ={'REQUEST_METHOD': 'POST',
                     'CONTENT_TYPE': self.headers['Content-Type'],
                     })

        query = QueryParams(dbnames=form["dbnames"].value.split(','),
                            dbkeys=form.getvalue("dbkeys", "").split(','),
                            version=int(form.getvalue("version", "1")),
                            smiles=form["smiles"].value.strip(),
                            return_count=int(form["return_count"].value),
                            similarity_cutoff=float(form["similarity_cutoff"].value)) #noqa

        if len(query.dbnames) != len(query.dbkeys):
            raise RuntimeError("Need key for each database.")

        return query

    def get_data(self, query, request_num, retry=0):
        global socket
        fp_binary, canon_smile = gpusim_utils.smiles_to_fingerprint_bin(
            query.smiles)
        fp_qba = QtCore.QByteArray(fp_binary)

        parameter_qba = QtCore.QByteArray()
        parameter_qds = QtCore.QDataStream(parameter_qba,
                                           QtCore.QIODevice.WriteOnly)

        parameter_qds.writeInt(len(query.dbnames))
        for name, key in zip(query.dbnames, query.dbkeys):
            parameter_qds.writeString(name.encode())
            parameter_qds.writeString(key.encode())

        parameter_qds.writeInt(request_num)
        parameter_qds.writeInt(query.return_count)
        parameter_qds.writeFloat(query.similarity_cutoff)
        parameter_qds << fp_qba

        socket.write(parameter_qba)
        socket.flush()

        response = QtCore.QSharedMemory(str(request_num))
        while not response.attach(QtCore.QSharedMemory.ReadOnly):
            if response.error() != QtCore.QSharedMemory.NotFound:
                # There are rare sporadic failures, likely due to a race
                # condition in QSharedMemory.  Bug has been observed online
                print(f"Request {request_num} failed", file=sys.stderr)
                print(response.error(), file=sys.stderr)
                if retry < MAX_RETRY:
                    return self.get_data(query, request_num, retry+1)
                else:
                    raise RuntimeError(response.errorString())

        returned_request = 0
        # Lock, check data, and unlock until it's populated
        start = time.time()
        while returned_request == 0:
            if not response.lock():
                raise RuntimeError(response.errorString())

            if (time.time() - start) > SERVER_RESULT_TIMEOUT:
                raise RuntimeError("Call timed out.")

            output_qba = QtCore.QByteArray(bytes(response.constData()))
            data_reader = QtCore.QDataStream(output_qba)
            returned_request = data_reader.readInt()
            response.unlock()

        response.detach()
        return output_qba

    def do_GET(self):
        self.send_error(404, 'Server unavailable.')

    def search_for_results(self, query):
        global search_mutex
        search_mutex.lock()
        try:
            request_num = random.randint(0, 2**31)
            output_qba = self.get_data(query, request_num)
            approximate_results, return_count, smiles, ids, scores = \
                self.deserialize_results(request_num, output_qba)

            return approximate_results, smiles, ids, scores
        finally:
            search_mutex.unlock()

    def flush_socket(self):
        global socket
        while not socket.atEnd():
            socket.readAll()

    def do_POST(self):

        if not self.path.startswith("/similarity_search_json"):
            return

        query = self.get_posted_query()

        if not check_socket():
            self.send_response(200)
            self.send_header('Content-type', 'text/json')
            self.end_headers()
            response = generate_error_response(query.version)
            self.wfile.write(json.dumps(response).encode('utf8'))
            return

        try:
            approx_results, smiles, ids, scores = \
                self.search_for_results(query)
        except RuntimeError as e:
            print(str(e), file=sys.stderr)  # Output to server log
            approx_results, smiles, ids, scores = SERVER_ERROR_RESULT

        self.do_json_POST(approx_results, smiles, ids, scores,
                          query.version)

    def deserialize_results(self, request_num, output_qba):
        data_reader = QtCore.QDataStream(output_qba)
        returned_request = data_reader.readInt()
        if request_num != returned_request:
            raise RuntimeError("Incorrect result ID returned!")
        return_count = data_reader.readInt()
        approximate_matches = data_reader.readUInt64()
        smiles, ids, scores = [], [], []
        for i in range(return_count):
            smiles.append(data_reader.readString().decode("utf-8"))
        for i in range(return_count):
            ids.append(data_reader.readString().decode("utf-8"))
        for i in range(return_count):
            scores.append(data_reader.readFloat())
        return approximate_matches, return_count, smiles, ids, scores

    def results2json(self, approx_results, smiles, ids, scores, version):
        results = None

        # Depending on the version of the API caller, we return results
        # in a different versioned format.
        if version == 1:
            results = []
            for cid, smi, score in zip(ids, smiles, scores):
                results.append([cid, smi, score])
        elif version == 2:
            results = {}
            results["approximate_count"] = approx_results
            results_list = []
            for cid, smi, score in zip(ids, smiles, scores):
                results_list.append([cid, smi, score])
            results["results"] = results_list

        return results

    def do_json_POST(self, approx_results, smiles, ids, scores, version):
        self.send_response(200)
        self.send_header('Content-type', 'text/json')
        self.end_headers()

        results_json = self.results2json(approx_results, smiles, ids, scores,
                                         version)
        self.wfile.write(json.dumps(results_json).encode('utf8'))


class GPUSimHTTPHandler(GPUSimHandler):
    _tmp_dir = tempfile.TemporaryDirectory()

    def _generateImageIfNeeded(self, fullpath):
        if (self.path.startswith('smiles_') and not os.path.exists(fullpath)):
            safe_smi = self.path.replace('smiles_', '').replace('.png', '')
            smi = safe_smi.replace('_-1-_', '/').replace(
                '_-2-_', '\\').replace('_-3-_', '#')
            gpusim_utils.smiles_to_image_file(smi, fullpath)

    def write_results_html(self, approx_results, src_smiles, smiles, scores,
                           ids):
        self.wfile.write(
            "Approximate Total Matching Compounds: {0}, returning {1}<p>".format( #noqa
                approx_results, len(smiles)).encode())
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

        query = self.get_posted_query()

        if not check_socket():
            self.send_response(200)
            self.send_header('Content-type', 'text/json')
            self.end_headers()
            response = generate_error_response(query.version)
            self.wfile.write(json.dumps(response).encode('utf8'))
            return

        try:
            approx_results, smiles, ids, scores = \
                self.search_for_results(query)
        except RuntimeError as e:
            print(str(e), file=sys.stderr)  # Output to server log
            approx_results, smiles, ids, scores = SERVER_ERROR_RESULT

        if self.path.startswith("/similarity_search_json"):
            self.do_json_POST(approx_results, smiles, ids, scores,
                              query.version)
        elif self.path.startswith("/similarity_search"):
            self.do_html_POST(approx_results, smiles, ids, scores,
                              query.smiles)

    def do_html_POST(self, approx_results, smiles, ids, scores, src_smiles):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

        f = open(os.path.join(SCRIPT_DIR, "index.html"), 'rb+')
        self.wfile.write(f.read())
        self.write_results_html(approx_results, src_smiles, smiles, scores,
                                ids)


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


def check_socket():
    global socket

    if socket is None:
        app = QtCore.QCoreApplication.instance()
        socket = QtNetwork.QLocalSocket(app)
    if not socket.isValid():
        socket_name = 'gpusimilarity'
        socket.connectToServer(socket_name)
    return socket.isValid()


def generate_error_response(version):
    error = None
    if version == 1:
        error = [['SERVER_REBOOTING', 'CC', '1.0']]
    elif version == 2:
        error = {"approximate_count": 0,
                 "error_code": "The server is currently restarting.",
                 "results": [['SERVER_REBOOTING', 'CC', '1.0']]}
    return error


def main():

    global server
    args = parse_args()

    # Create the QApplication we run as
    QtCore.QCoreApplication([])

    # Start the GPU backend
    cmdline = [GPUSIM_EXEC]
    if args.cpu_only:
        cmdline.append('--cpu_only')
    cmdline += ['--gpu_bitcount', args.gpu_bitcount]
    cmdline += args.dbnames
    backend_proc = subprocess.Popen(cmdline)

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
