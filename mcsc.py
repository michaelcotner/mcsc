VERSION = '0.1'

import urllib.request

def _get_version():
    # hosted file to compare to
    url = 'https://raw.githubusercontent.com/michaelcotner/mcsc/main/mcsc.py'

    try:
        with urllib.request.urlopen(url) as response:
            # get version from the first line
            file_head = response.readline().decode('utf-8')
            if file_head.startswith('VERSION = '):
                version = file_head[len('VERSION = '):].strip('\'')
                return version
    except Exception as e:
        print('Failed to get version:', e)
        return None
    
    print('Failed to get version')
    return None


def test1():
    return 'test1'

def test():
    return 'test2'