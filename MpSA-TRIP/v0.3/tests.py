#!/usr/bin/env python
# encoding: utf-8

import trip_0_3 as trip
import unittest, os

class TestingLoadMainConfig(unittest.TestCase):

    def test_true_config(self):
        config = trip.load_main_config("./test/config.json")
        self.assertEqual(type(config), dict)

    def test_length_config(self):
        config = trip.load_main_config("./test/config.json")
        self.assertEqual(len(config), 6)

    def test_false_config(self):
        self.assertRaises(IOError, trip.load_main_config, 'config.json')

class TestingLoadLoggingConfig(unittest.TestCase):

    def test_true_path(self):
        self.assertTrue(os.path.exists("logging.json"))

# class TestOpenInput(unittest.TestCase):

#     def test_true_input(self):
#         INPUT = "testing_file.fastq.gz"


if __name__ == '__main__':
    unittest.main()