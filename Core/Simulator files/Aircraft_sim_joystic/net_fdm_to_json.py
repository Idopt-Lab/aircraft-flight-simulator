"""The net_fdm to json module

This module accepts incoming udp scoket traffic and forwards it to a websocket
"""

import struct


def listen_matlab(format, fields, dims, sock_incoming):
    """Listen to Matlab UDP traffic and convert to json"""
    try:
        data, _ = sock_incoming.recvfrom(8192)
        unpacked_data = struct.unpack(format, data)
        json_data = {}
        idy = 0
        for idx, field in enumerate(fields):
            json_data[field] = list(unpacked_data[idy:idy+dims[idx]])
            idy += dims[idx]
        return json_data
    except Exception:
        #print("Problem recieving from Matlab.")
        return None