import lxml.etree as lxml
import argparse 
import sys

def Parser(optionfile,entry):
    parser=lxml.XMLParser(remove_comments=True)
    tree = lxml.parse(optionfile,parser)
    root = tree.getroot() 
        
    return root.find(entry)

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
