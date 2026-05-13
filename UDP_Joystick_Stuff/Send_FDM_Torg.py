# Send_FDM.py
# HOTAS Warthog stick + Warthog throttle via pygame
# T-Rudder pedals via /dev/input/js1 (Linux joystick API) to bypass SDL device listing issues.

import os
os.environ["SDL_JOYSTICK_ALLOW_BACKGROUND_EVENTS"] = "1"

import socket
import struct
import pygame


# =========================
# UDP CONFIG
# =========================
UDP_IP = "127.0.0.1"
UDP_PORT = 5505
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)


# =========================
# UI PRINT HELPER
# =========================
class TextPrint:
    def __init__(self):
        self.reset()
        self.font = pygame.font.Font(None, 25)

    def tprint(self, screen, text):
        text_bitmap = self.font.render(text, True, (0, 0, 0))
        screen.blit(text_bitmap, (self.x, self.y))
        self.y += self.line_height

    def reset(self):
        self.x = 10
        self.y = 10
        self.line_height = 15


def apply_deadzone(x, dz=0.05):
    return 0.0 if abs(x) < dz else x


def norm_name(s: str) -> str:
    return (s or "").lower().replace("-", " ").replace("_", " ").strip()


def dump_devices(joysticks):
    print("\n=== DETECTED DEVICES (pygame) ===")
    for jid, joy in joysticks.items():
        print(
            f"[id={jid}] name='{joy.get_name()}' "
            f"axes={joy.get_numaxes()} buttons={joy.get_numbuttons()} hats={joy.get_numhats()}"
        )
    print("=================================\n")


def find_best_device(joysticks, include_keywords, exclude_keywords=None):
    exclude_keywords = exclude_keywords or []
    best = None
    best_score = -1

    for joy in joysticks.values():
        name = norm_name(joy.get_name())
        if any(k in name for k in exclude_keywords):
            continue
        score = sum(1 for k in include_keywords if k in name)
        if score > best_score and score > 0:
            best = joy
            best_score = score

    return best


def choose_throttle(throttle_dev, deadzone):
    """
    Prefer axes 2 & 3 (avg). Fallback to axes 0 & 1 (avg). Then axis 0.
    Returns (throttle_value, mode_string).
    """
    if throttle_dev is None:
        return 0.0, "none"

    na = throttle_dev.get_numaxes()

    def dz(v): return apply_deadzone(v, deadzone)

    if na >= 4:
        a2 = dz(throttle_dev.get_axis(2))
        a3 = dz(throttle_dev.get_axis(3))
        return 0.5 * (a2 + a3), "avg(thr axes 2,3)"

    if na >= 2:
        a0 = dz(throttle_dev.get_axis(0))
        a1 = dz(throttle_dev.get_axis(1))
        return 0.5 * (a0 + a1), "avg(thr axes 0,1)"

    if na >= 1:
        return dz(throttle_dev.get_axis(0)), "thr axis 0"

    return 0.0, "none"


class LinuxJoystickReader:
    """
    Non-blocking reader for /dev/input/jsX (Linux joystick API).
    Reads axis values like jstest does.
    """
    JS_EVENT_AXIS = 0x02
    JS_EVENT_INIT = 0x80

    def __init__(self, path: str):
        self.path = path
        self.fd = os.open(path, os.O_RDONLY | os.O_NONBLOCK)
        self.axis = {}  # axis_index -> value in [-1, 1]

    def close(self):
        try:
            os.close(self.fd)
        except Exception:
            pass

    def poll(self):
        fmt = "IhBB"  # time(uint32), value(int16), type(uint8), number(uint8)
        size = struct.calcsize(fmt)

        while True:
            try:
                data = os.read(self.fd, size)
                if len(data) < size:
                    break

                _, value, etype, number = struct.unpack(fmt, data)
                etype = etype & ~self.JS_EVENT_INIT

                if etype == self.JS_EVENT_AXIS:
                    self.axis[number] = max(-1.0, min(1.0, value / 32767.0))

            except BlockingIOError:
                break
            except OSError:
                break

    def get_axis(self, idx: int, default: float = 0.0) -> float:
        return float(self.axis.get(idx, default))


def main():
    pygame.init()
    pygame.joystick.init()

    screen = pygame.display.set_mode((560, 760))
    pygame.display.set_caption("HOTAS -> UDP FDM (Stick+Throttle pygame, Pedals js1)")

    clock = pygame.time.Clock()
    text_print = TextPrint()

    joysticks = {}

    DEADZONE = 0.05
    FPS = 60

    # Debug axis printing
    last_axes = {}
    AXIS_PRINT_THRESHOLD = 0.20

    # Pedals via Linux js API (your T-Rudder is /dev/input/js1)
    pedals_path = "/dev/input/js1"
    pedals_js = None
    pedals_ok = False
    try:
        pedals_js = LinuxJoystickReader(pedals_path)
        pedals_ok = True
        print(f"Pedals reader opened: {pedals_path}")
    except Exception as e:
        print(f"WARNING: Could not open {pedals_path}. Yaw will be 0. Error: {e}")

    pygame.event.pump()

    done = False
    while not done:
        # --- events ---
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                done = True

            if event.type == pygame.JOYDEVICEADDED:
                joy = pygame.joystick.Joystick(event.device_index)
                joy.init()
                joysticks[joy.get_instance_id()] = joy
                print(
                    f"Joystick connected: id={joy.get_instance_id()} name='{joy.get_name()}' "
                    f"axes={joy.get_numaxes()} buttons={joy.get_numbuttons()}"
                )
                dump_devices(joysticks)

            if event.type == pygame.JOYDEVICEREMOVED:
                if event.instance_id in joysticks:
                    name = joysticks[event.instance_id].get_name()
                    del joysticks[event.instance_id]
                    print(f"Joystick disconnected: id={event.instance_id} name='{name}'")
                    dump_devices(joysticks)

            if event.type == pygame.JOYBUTTONDOWN:
                print(f"Joystick button {event.button} pressed (id={event.instance_id})")

        # --- choose devices ---
        stick = find_best_device(
            joysticks,
            include_keywords=["warthog", "joystick", "thustmaster", "thrustmaster"],
            exclude_keywords=["throttle", "rudder", "pedal"]
        )

        throttle_dev = find_best_device(
            joysticks,
            include_keywords=["warthog", "throttle"],
            exclude_keywords=["rudder", "pedal", "joystick"]
        )

        # --- read axes ---
        roll = 0.0
        pitch = 0.0
        yaw = 0.0
        throttle = 0.0
        thr_mode = "none"

        if stick and stick.get_numaxes() >= 2:
            roll = apply_deadzone(stick.get_axis(0), DEADZONE)
            pitch = apply_deadzone(stick.get_axis(1), DEADZONE)

        throttle, thr_mode = choose_throttle(throttle_dev, DEADZONE)

        # Flip throttle sign if you want positive-up:
        throttle = -throttle

        if pedals_ok and pedals_js is not None:
            pedals_js.poll()
            yaw = apply_deadzone(pedals_js.get_axis(2, 0.0), DEADZONE)  # T-Rudder yaw = axis 2

        # --- debug axis prints ---
        for jid, joy in joysticks.items():
            for ai in range(joy.get_numaxes()):
                v = joy.get_axis(ai)
                key = (jid, ai)
                v0 = last_axes.get(key, v)
                if abs(v - v0) > AXIS_PRINT_THRESHOLD:
                    print(f"AXIS MOVE: id={jid} name='{joy.get_name()}' axis={ai} value={v:+.3f}")
                last_axes[key] = v

        # --- UDP send ---
        msg = struct.pack("ffff", float(roll), float(pitch), float(yaw), float(throttle))
        sock.sendto(msg, (UDP_IP, UDP_PORT))

        # --- UI ---
        screen.fill((255, 255, 255))
        text_print.reset()

        text_print.tprint(screen, f"UDP Target: {UDP_IP}:{UDP_PORT}")
        text_print.tprint(screen, f"Pygame Devices: {len(joysticks)}")
        text_print.tprint(screen, "")
        text_print.tprint(screen, f"Stick:    {(stick.get_name() if stick else 'None')}")
        text_print.tprint(screen, f"Throttle: {(throttle_dev.get_name() if throttle_dev else 'None')}")
        text_print.tprint(screen, f"Pedals:   {'/dev/input/js1' if pedals_ok else 'None'}")
        text_print.tprint(screen, "")
        text_print.tprint(screen, f"Roll:     {roll:+.3f}  (stick axis 0)")
        text_print.tprint(screen, f"Pitch:    {pitch:+.3f}  (stick axis 1)")
        text_print.tprint(screen, f"Yaw:      {yaw:+.3f}  (js1 axis 2)")
        text_print.tprint(screen, f"Throttle: {throttle:+.3f}  ({thr_mode}, inverted)")
        text_print.tprint(screen, "")

        if not pedals_ok:
            text_print.tprint(screen, "Could not open /dev/input/js1 for pedals.")
            text_print.tprint(screen, "Check: ls -l /dev/input/js1 (permissions)")

        pygame.display.flip()
        clock.tick(FPS)

    if pedals_js is not None:
        pedals_js.close()
    pygame.quit()


if __name__ == "__main__":
    main()

