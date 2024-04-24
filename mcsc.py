VERSION = '0.1'

import urllib.request

def _get_version_from_github():
    # URL of the raw Python file on GitHub
    github_url = 'https://raw.githubusercontent.com/michaelcotner/mcsc/main/mcsc.py'

    # Fetch the raw content of the Python file
    try:
        with urllib.request.urlopen(github_url) as response:
            # Extract the value of the VERSION variable from the first line
            file_head = response.readline().decode('utf-8')
            if file_head.startswith('VERSION = '):
                version = file_head[len('VERSION = '):].strip()
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