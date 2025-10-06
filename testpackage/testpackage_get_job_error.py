
import os

from argparse import ArgumentParser


parser = ArgumentParser(prog="Log file parser"
                        ,description="Parses a log file and extracts job error messages."
                        )

parser.add_argument("log_file",help="Log file to parse")

args= parser.parse_args()

file= args.log_file

with open(file,'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.split(" ")[0]=="EXIT_CODE_FROM_JOB":
            raise SystemExit(int(line.split(" ")[1].rstrip("\n")))
