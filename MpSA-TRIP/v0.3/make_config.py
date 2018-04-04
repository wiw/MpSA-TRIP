#!/usr/bin/env python
# encoding: utf-8

"""
    'make_config.py'
    Interactive config generator from data avalable from user input.
    Used to automatically create a configuration file for the program 'trip_0.3.py'
    author: Anton V. Ivankin
    e-mail: anton.ivankin@gmail.com
    source: https://github.com/wiw/MpSA-TRIP/
"""

import os, json, logging

logging.basicConfig(filename=None, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%Y.%m.%d %H:%M:%S')
Logger = logging.getLogger(__name__)

path_to_config = os.path.join(os.getcwd(), "config.json")

def save_config(user_input, path):
    try:
        with open(path, "wb") as handle:
            json.dump(user_input, handle)
            Logger.info("Config are writing succesfully to '{}'".format(path))
    except:
        Logger.exeption("Doesn't save config file by path '{}'".format(path))

def user_input_assembly():
    pass

def main():
    config_dict = user_input_assembly()
    save_config(config_dict, path_to_config)

if __name__ == '__main__':
    main()