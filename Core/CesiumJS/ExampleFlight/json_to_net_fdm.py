"""The json to net_fdm module

This module accepts incoming websocket traffic and forwards it to a udp socket
"""


import struct


def send_matlab(json_data, format, fields, dims, sock_outgoing, MATLAB):
    """Convert json to net_fdm and send to Matlab"""
    #json_data = websocket_data
    values = []
    for idx, field in enumerate(fields):
        vector = json_data.get(field, [0.0] * dims[idx])
        values.extend(vector) # Convert the json into an array of 
    net_fdm = struct.pack(format, *values)
    sock_outgoing.sendto(net_fdm, MATLAB)