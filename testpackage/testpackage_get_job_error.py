
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
        error_code = line.split(" ")
        if error_code[0]=="EXIT_CODE_FROM_JOB" and error_code[1].rstrip("\n")!="0":
            raise SystemExit(int(line.split(" ")[1].rstrip("\n")))


