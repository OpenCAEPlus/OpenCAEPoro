# coding=utf8
# the above tag defines encoding for this document and is for Python 2.x compatibility

import re
import os
import sys

pathname = sys.argv[1]
print(pathname)

regex1 = r"(?<=[\n\r])^#"
regex2 = r"^-{4,}.+$(?![\r\n])"
regex3 = r"^"
header = """---
head:
    - - link
      - rel: stylesheet
        href: https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/katex.min.css
---

"""

dir_list = os.listdir(pathname)

for file in dir_list:
    if file.endswith(".md"):
        print(file)
        filename = os.path.join(pathname,file)
        f = open(filename,"r")
        test_str = f.read()
        f.close()
        result = re.sub(regex1, "##", test_str, 0, re.MULTILINE)
        result = re.sub(regex2, "", result, 0, re.DOTALL | re.MULTILINE)
        result = re.sub(regex3, header, result, 1)

        f = open(filename,"w")
        f.write(result)
        f.close()
        