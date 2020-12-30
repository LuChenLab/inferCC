#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
2018.12.25

Move data from this computer to disk
"""
import argparse as ap
import os
import sys
import paramiko


def iter_all_over(target, output, sftp):
    u"""
    iterover all directories and files to get the path of html file
    """
    data = [str(x) for x in sftp.listdir_iter(target)]

    for i in data:
        file_name = i.split()[-1]
        if "d" in i.split()[0]:
            iter_all_over(os.path.join(target, file_name), output, sftp)
        else:
            if target.endswith("outs") and file_name.endswith("web_summary.html"):
                file_prefix = os.path.basename(os.path.dirname(target))

                print(file_prefix, file_name)

                sftp.get(
                    remotepath=os.path.join(target, file_name),
                    localpath=os.path.join(output, file_prefix + "_" + file_name)
                )


def main(args):
    u"""
    Main functions
    :param args: output
    :return:
    """
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    iter_all_over(args.dir, args.output, sftp)


if __name__ == '__main__':

    parser = ap.ArgumentParser("Pull down the web_summary.html from server")

    parser.add_argument(
        "-host",
        help="Host name",
        type=str,
        required=True
    )

    parser.add_argument(
        "-port",
        help="Port",
        type=int,
        required=True,
        default=22
    )

    parser.add_argument(
        "-user",
        help="User name",
        type=str,
        required=True
    )

    parser.add_argument(
        "-passwd",
        help="Password",
        type=str,
        required=True
    )

    parser.add_argument(
        "-dir",
        help="Target directory",
        type=str,
        required=True
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        help="Path to output directory"
    )


    try:
        if len(sys.argv) <= 1:
            parser.print_help()
        else:
            args = parser.parse_args(sys.argv[1:])

            t = paramiko.Transport((args.host, args.port))

            t.connect(username=args.user, password=args.passwd)
            sftp = paramiko.SFTPClient.from_transport(t)

            remote = args.remote

            main(args)
    except Exception as err:
        print(err.__traceback__)
    finally:
        sftp.close()
