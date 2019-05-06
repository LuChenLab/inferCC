#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@Created at 2019.05.06

This scripts is used to pull down all the loom files from server
"""
import os
import sys

import logging

import json

from argparse import ArgumentParser, ArgumentError, Namespace
import paramiko


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__dir__ = os.path.abspath(os.path.dirname(__file__))

__author__ = "ygidtu"
__email__ = "ygidtu@gmail.com"
__version__ = "0.1"


remote = "/mnt/raid62/Lung_cancer_10x/CellRanger/lung_cancer"


def connect(host: str, port: int, user: str, passwd: str) -> paramiko.SFTP:
    t = paramiko.Transport((host, port))

    t.connect(username=user, password=passwd)
    sftp = paramiko.SFTPClient.from_transport(t)

    return sftp


def iter_all_over(sftp: paramiko.SFTP, target: str, postfix: str) -> list:
    u"""
    iter over all directories and files to get the path of html file
    :param sftp:
    :param target: path to target directory
    :param postfix: postfix of target files
    :return list of path
    """
    results = []
    data = [str(x) for x in sftp.listdir_iter(target)]
    for i in data:
        file_name = i.split()[-1]
        if "d" in i.split()[0]:
            results += iter_all_over(
                sftp=sftp,
                target=os.path.join(target, file_name),
                postfix=postfix
            )
        else:
            if postfix != "*" and i.endswith(postfix):
                results.append(os.path.join(target, file_name))
    return results


def main(args: Namespace):
    u"""
    Main function as the entry point of this scripts
    :param args: parsed argparse
    :return:
    """
    logger.info(args)
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    temp_progress = os.path.join(args.out, "progress.json")

    logger.info("Connect to {0}:{1}".format(args.host, args.port))
    sftp = connect(
        host=args.host,
        port=args.port,
        user=args.user,
        passwd=args.passwd
    )

    if not args.force and os.path.exists(temp_progress):
        with open(temp_progress) as r:
            files = json.load(r)
    else:
        logger.info("Finding targets")
        files = iter_all_over(
            sftp=sftp,
            target=args.dir,
            postfix=args.file
        )

        with open(temp_progress, "w+") as w:
            json.dump(files, w, indent=4)

    logger.info("Pulling")
    for i in files:

        local_path = os.path.join(args.out, i.replace(args.dir, "").lstrip("/"))
        logger.info("Getting {0} to {1}".format(i, local_path))

        if not os.path.exists(os.path.dirname(local_path)):
            os.makedirs(os.path.dirname(local_path))

        sftp.get(
            remotepath=i,
            localpath=local_path
        )


if __name__ == '__main__':
    parser = ArgumentParser(description="Pull all loom files from server")

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__
    )

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
        "-out",
        help="Output directory",
        type=str,
        default=__dir__
    )

    parser.add_argument(
        "-file",
        help="postfix of target files",
        type=str,
        default="*"
    )

    parser.add_argument(
        "-force",
        help="Force to get files",
        action="store_true"
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        exit(0)

    try:
        args = parser.parse_args(sys.argv[1:])
        main(args)
    except ArgumentError as err:
        print(err)
        parser.print_help()

