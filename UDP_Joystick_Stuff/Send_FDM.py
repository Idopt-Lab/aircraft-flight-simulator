# Send FDM
# For validating sending joystick data along FDM port.
# Written by Casey Chamberlain
# 11-13-2025

# TURN OFF DEADZONE FOR ZERO CONTROL INPUT 
# Do a mock controller like Ben has
# Post stuff to github

# Utilities
import os

os.environ["SDL_JOYSTICK_ALLOW_BACKGROUND_EVENTS"]="1"

import pygame
pygame.init()
import socket
import struct
import time
from pygame import *

# Initialize UDP port
UDP_IP = "127.0.0.1"
UDP_PORT = 5505 # Is this supposed to be the same as the receiver?
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM) # PROGRAMMING SOCKS HECK YEAH

# Code begins here
# Source: https://www.pygame.org/docs/ref/joystick.html#pygame.joystick.Joystick.init
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

    def indent(self):
        self.x += 10

    def unindent(self):
        self.x -= 10


def main():
    # Set the width and height of the screen (width, height), and name the window.
    screen = pygame.display.set_mode((500, 700))
    pygame.display.set_caption("Joystick example")

    # Used to manage how fast the screen updates.
    clock = pygame.time.Clock()

    # Get ready to print.
    text_print = TextPrint()

    # pygame.event.set_grab()

    # This dict can be left as-is, since pygame will generate a
    # pygame.JOYDEVICEADDED event for every joystick connected
    # at the start of the program.
    joysticks = {}

    # pygame.event.set_grab(False)

    done = False
    # Variables
    pitch = 0
    roll = 0
    yaw = 0
    throttle = 0
    while not done:
        # Event processing step.
        # Possible joystick events: JOYAXISMOTION, JOYBALLMOTION, JOYBUTTONDOWN,
        # JOYBUTTONUP, JOYHATMOTION, JOYDEVICEADDED, JOYDEVICEREMOVED
        
        for event in pygame.event.get(): # This section runs ONLY if there's input!
            if event.type == pygame.QUIT:
                done = True  # Flag that we are done so we exit this loop.

            if event.type == pygame.JOYBUTTONDOWN: # Modify this to identify which button is pressed
                print(f"Joystick button {event.button} pressed.")
                if event.button == 0:
                    joystick = joysticks[event.instance_id]
                    if joystick.rumble(0, 0.7, 500):
                        print(f"Rumble effect played on joystick {event.instance_id}")

            if event.type == pygame.JOYBUTTONUP:
                print(f"Joystick button {event.button} released.")

            # Handle hotplugging
            if event.type == pygame.JOYDEVICEADDED:
                # This event will be generated when the program starts for every
                # joystick, filling up the list without needing to create them manually.
                joy = pygame.joystick.Joystick(event.device_index)
                joysticks[joy.get_instance_id()] = joy
                print(f"Joystick {joy.get_instance_id()} connencted")

            if event.type == pygame.JOYDEVICEREMOVED:
                del joysticks[event.instance_id]
                print(f"Joystick {event.instance_id} disconnected")

            if event.type == pygame.JOYAXISMOTION:
                # Notes:
                # joystick axis 0, 1 (left stick x y)
                # joystick axis 2, 3 (right stick x y)
                roll = joystick.get_axis(2)
                pitch = joystick.get_axis(3)
                yaw = joystick.get_axis(0)
                throttle = joystick.get_axis(1)
                msg = struct.pack('ffff', roll, pitch, yaw, throttle)
                print(f"Pitch input: {pitch}")
                print(f"Roll input: {roll}")
                print(f"Yaw input: {yaw}")
                print(f"Throttle input: {throttle}")
                sock.sendto(msg, (UDP_IP, UDP_PORT))

            # Adding support for MnK
            # if event.type == pygame.KEYDOWN:
            #     # ROLL
            #     if event.key == K_a:
            #         print("K_a")
            #         roll = -1
            #     elif event.key == K_d:
            #         print("K_d")
            #         roll = 1
            #     else:
            #         roll = 0
            #     # PITCH
            #     if event.key == K_w:
            #         print("K_w")
            #         pitch = -1
            #     elif event.key == K_s:
            #         print("K_s")
            #         pitch = 1
            #     else:
            #         pitch = 0
            #     # YAW
            #     if event.key == K_e:
            #         print("K_e")
            #         yaw = 1
            #     elif event.key == K_q:
            #         print("K_q")
            #         yaw = -1
            #     else:
            #         yaw = 0
            #     # THROTTLE
            #     if event.key == K_LSHIFT:
            #         print("K_lshift")
            #         throttle = throttle + 0.1
            #     elif event.key == K_LCTRL:
            #         print("K_lctrl")
            #         throttle = throttle - 0.1
            #     else:
            #         throttle = throttle
            pressed_key = pygame.key.get_pressed()
            released_key = pygame.key.get_just_released() 
            # ROLL
            if pressed_key[pygame.K_e]:
                print("Key e is being pressed.")
                roll = 1 # pos is right/cw
            if pressed_key[pygame.K_q]:
                print("Key q is being pressed.")
                roll = -1 # neg is left/ccw
            if released_key[pygame.K_e] | released_key[pygame.K_q]:
                print("Roll key has been released.")
                roll = 0
            # PITCH
            if pressed_key[pygame.K_w]:
                print("Key w is being pressed.")
                pitch = -1
            if pressed_key[pygame.K_s]:
                print("Key s is being pressed.")
                pitch = 1
            if released_key[pygame.K_w] | released_key[pygame.K_s]:
                print("Pitch key has been released.")
                pitch = 0
            # YAW
            if pressed_key[pygame.K_d]:
                print("Key d is being pressed.")
                yaw = 1
            if pressed_key[pygame.K_a]:
                print("Key a is being pressed.")
                yaw = -1
            if released_key[pygame.K_a] | released_key[pygame.K_d]:
                print("Yaw keys released.")
                yaw = 0
            # THROTTLE
            if pressed_key[pygame.K_LSHIFT]:
                print("Increasing throttle.")
                throttle = throttle + 0.5
            if pressed_key[pygame.K_LCTRL]:
                print("Decreasing throttle.")
                throttle = throttle - 0.5
            if released_key[pygame.K_LSHIFT] | released_key[pygame.K_LCTRL]:
                print("Throttle released.")
                throttle = throttle

            # msg = struct.pack('ffff', roll, pitch, yaw, throttle)
            # print(f"Roll input: {roll}")
            # sock.sendto(msg, (UDP_IP, UDP_PORT))

            # Exit thing
            if pressed_key[pygame.K_ESCAPE]:
                pygame.quit()

        # Drawing step
        # First, clear the screen to white. Don't put other drawing commands
        # above this, or they will be erased with this command.
        screen.fill((255, 255, 255))
        text_print.reset()

        # Get count of joysticks.
        joystick_count = pygame.joystick.get_count()

        text_print.tprint(screen, f"Number of joysticks: {joystick_count}")
        text_print.indent()

        # if joystick_count < 1, use MnK
        if joystick_count < 1:
            text_print.tprint(screen, f"Using Mouse and Keyboard.")

            # Print each key input
            text_print.tprint(screen, f"Up: {K_w}")

        else:

            # For each joystick:
            for joystick in joysticks.values():
                jid = joystick.get_instance_id()

                text_print.tprint(screen, f"Joystick {jid}")
                text_print.indent()

                # Get the name from the OS for the controller/joystick.
                name = joystick.get_name()
                text_print.tprint(screen, f"Joystick name: {name}")

                guid = joystick.get_guid()
                text_print.tprint(screen, f"GUID: {guid}")

                power_level = joystick.get_power_level()
                text_print.tprint(screen, f"Joystick's power level: {power_level}")

                # Usually axis run in pairs, up/down for one, and left/right for
                # the other. Triggers count as axes.
                axes = joystick.get_numaxes()
                text_print.tprint(screen, f"Number of axes: {axes}")
                text_print.indent()

                for i in range(axes):
                    axis = joystick.get_axis(i)
                    text_print.tprint(screen, f"Axis {i} value: {axis:>6.3f}")
                text_print.unindent()

                buttons = joystick.get_numbuttons()
                text_print.tprint(screen, f"Number of buttons: {buttons}")
                text_print.indent()

                for i in range(buttons):
                    button = joystick.get_button(i)
                    text_print.tprint(screen, f"Button {i:>2} value: {button}")
                text_print.unindent()

                hats = joystick.get_numhats()
                text_print.tprint(screen, f"Number of hats: {hats}")
                text_print.indent()

                # Hat position. All or nothing for direction, not a float like
                # get_axis(). Position is a tuple of int values (x, y).
                for i in range(hats):
                    hat = joystick.get_hat(i)
                    text_print.tprint(screen, f"Hat {i} value: {str(hat)}")
                text_print.unindent()

                text_print.unindent()

        # Go ahead and update the screen with what we've drawn.
        pygame.display.flip()

        # Limit to 30 frames per second.
        clock.tick(30)


if __name__ == "__main__":
    main()
    # If you forget this line, the program will 'hang'
    # on exit if running from IDLE.
    pygame.quit()