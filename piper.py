#### This little script filters the CIGAR strings and only when the alignment has a mimumum amount of bases it will be pushed to stdout
import sys
import re
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

import optparse

parser = optparse.OptionParser()

parser.add_option('-m', '--minAlignment',
    action="store", dest="minAlignment",
    help="set the minimal alignment length", default="spam")

options, args = parser.parse_args()


if __name__ == "__main__":
    for line in sys.stdin:
		if line[0:1] == "@":
			sys.stdout.write(line)
			continue
#		sys.stderr.write("DEBUG: got line: " + line)
        #sys.stdout.write(line)
        



		regex = r"(\d+)M"

		test_str = line.split('\t')[5]

		matches = re.finditer(regex, test_str)

		alignMatch = 0
		for matchNum, match in enumerate(matches):
		    matchNum = matchNum + 1
		    
#		    print ("Match {matchNum} was found at {start}-{end}: {match}".format(matchNum = matchNum, start = match.start(), end = match.end(), match = match.group()))
		    
		    for groupNum in range(0, len(match.groups())):
		        groupNum = groupNum + 1
		        
#		        print ("Group {groupNum} found at {start}-{end}: {group}".format(groupNum = groupNum, start = match.start(groupNum), end = match.end(groupNum), group = match.group(groupNum)))
		        alignMatch += int(match.group(groupNum))

		if alignMatch > int(options.minAlignment):
			sys.stdout.write(line)
#		print alignMatch

