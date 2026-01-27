"""Main Python module

Runs X3D example with Matlab integration
"""

import threading
import queue
import json
import socket
from json_to_net_fdm import send_matlab
from net_fdm_to_json import listen_matlab
from http.server import HTTPServer, SimpleHTTPRequestHandler
from websockets.sync.server import serve


WEBSERVER=('localhost', 8000)
WEBSOCKET = ('localhost', 8765)
MATLAB_INCOMING = ('localhost', 50000)
MATLAB_OUTGOING = ('localhost', 50001)

exit_event = threading.Event() # Press enter to safely end session
ws_incoming = queue.Queue()  # Messages received from websocket clients
ws_outgoing = queue.Queue()  # Messages to send to websocket clients
sock_incoming = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock_outgoing = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

def create_packet_formats(filename):
    """Take a ...format.txt and convert into a usable format string and json format"""
    fields = [] # json fields
    dims = [] # json dimensions
    format = '>' # little endian
    with open(filename, 'r') as f: # https://sourceforge.net/p/flightgear/flightgear/ci/next/tree/src/Network/net_fdm.hxx
        for line in f:
            name, type_str = line.split()
            format += type_str # create struct format string
            fields.append(name) # create array of json fields
            dims.append(len(type_str)) # create array of vector sizes
            #fields[name] = [-1.0] * dim # create empty json container
    return format, fields, dims

def run_http_server(format, fields, dims):
    """Run the webserver for the X3D"""
    httpd = HTTPServer(WEBSERVER, SimpleHTTPRequestHandler)
    httpd.timeout = 1
    print(f"Running Webserver at http://{WEBSERVER[0]}:{WEBSERVER[1]}")
    while not exit_event.is_set():
        httpd.handle_request()
    #print("close http loop")

def run_main_loop(format, fields, dims):
    """A syncronous implementation of 2-way communication bridge"""
    
    sock_incoming.bind(MATLAB_INCOMING)
    sock_incoming.setblocking(False)
    while not exit_event.is_set():
        matlab_data = listen_matlab(format, fields, dims, sock_incoming)
        if matlab_data is not None: # If matlab sent data then send a websocket
            ws_outgoing.put(matlab_data)

        try: # If websocket sent data then send it to matlab
            websocket_data = ws_incoming.get_nowait()
            send_matlab(websocket_data, format, fields, dims, sock_outgoing, MATLAB_OUTGOING)
        except queue.Empty:
            pass
    #print("close main loop")

def handler(websocket):
    while not exit_event.is_set():
        #try:
        #    message = websocket.recv(timeout=0.1)
        #    if message:
        #        ws_incoming.put(json.loads(message))
        #except TimeoutError:
        #    pass

        try:
            outgoing_data = ws_outgoing.get_nowait()
            #print("Type being sent:", type(outgoing_data))
            websocket.send(json.dumps(outgoing_data))
            #print(json.dumps(outgoing_data))
        except queue.Empty:
            pass
    #print("close websocket server")

def run_websocket_server(websocket_server):
    with websocket_server as Server:
        Server.serve_forever()

if __name__ =="__main__":
    websocket_server = serve(handler, WEBSOCKET[0], WEBSOCKET[1])

    format, fields, dims = create_packet_formats("net_fdm_format.txt")
    t1 = threading.Thread(target=run_http_server, args=(format, fields, dims,), daemon=True)
    t2 = threading.Thread(target=run_main_loop, args=(format, fields, dims,), daemon=True)
    t3 = threading.Thread(target=run_websocket_server, args=(websocket_server,), daemon=True)

    t1.start()
    t2.start()
    t3.start()

    input("Press [Enter] to stop the program...\n")
    print("Shutting down...")
    exit_event.set()
    websocket_server.shutdown()

    # Give threads time to clean up
    t1.join()
    t2.join()
    t3.join()
    print("Threads closed.")