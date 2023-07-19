import sys
import re

arg = sys.argv[1]

matches = re.findall(r'(\d+)', arg)
for match in matches:
    arg = arg.replace(match, str(int(match) - 747))

print(arg)
