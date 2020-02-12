#!/usr/bin/env python3
import argparse
import sys
sys.path.insert(0, '~/Desktop/codes/')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = ' ')
    parser.add_argument('-a', required= True, help=' ')
    parser.add_argument('-p', default='none', nargs='*', required=True, help=' ')

    ns = parser.parse_args()



