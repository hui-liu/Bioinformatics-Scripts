import sys
import os

def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0

def file_size(file_path):
    """
    this function will return the file size
    """
    if os.path.isfile(file_path):
        file_info = os.stat(file_path)
        return convert_bytes(file_info.st_size)

for file_path in os.listdir(sys.argv[1]):
    if file_path[-2:] == "gz": continue
    try:
        s, u = file_size(sys.argv[1] + "/" + file_path).split()
    except AttributeError:
        print "skip the dir %s" % file_path
    s = float(s)
    if u in ['GB', 'TB']:
        #print file_path, s, u
        os.system('gzip %s' % sys.argv[1] + "/" + file_path)
    elif u == "MB" and s > 50:
        #print file_path, s, u
        os.system('gzip %s' % sys.argv[1] + "/" + file_path)
