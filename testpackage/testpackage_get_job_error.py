

from argparse import ArgumentParser


parser = ArgumentParser(prog="Log file parser"
                        ,description="Parses a log file and extracts job error messages."
                        )

parser.add_argument("log_file",help="Log file to parse")

args= parser.parse_args()

file= args.log_file

hostname_count=0
hostname_target=-1


with open(file,'r') as f:
    lines = f.readlines()
    for line in lines:
        line_list = line.split(" ")
        if line_list[0]=="EXIT_CODE_FROM_JOB" and line_list[1].rstrip("\n")!="0":
            raise SystemExit(int(line_list[1].rstrip("\n")))
        elif line_list[0]=="HOSTNAME":
            hostname_target=int(line_list[2].rstrip("\n"))
            hostname_count+=1
    f.close()

if hostname_target!=hostname_count:
    raise SystemExit(f"::error::{hostname_target-hostname_count} node(s) likely failed silently")

