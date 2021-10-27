from gooey import Gooey, GooeyParser
from version import __version__ as version
import sys

def ArgParse(defaults):
    """
    Asks user for all variables integrated into a GUI.
    """
    # Asks user for all variables
    # Graphical user interface
    
    parser = GooeyParser()
    
    # Analysis
    exp_group = parser.add_argument_group('Setup')
    exp_group.add_argument("xname", metavar= 'Experiment name', type = str, default= defaults['xname'])
    exp_group.add_argument("outdir", metavar = 'Output directory', type = str, widget="DirChooser")
    exp_group.add_argument('rawdata', metavar ='Raw data file', type = str, widget = 'FileChooser', help='.csv exported from the RT-PCR Software')
    exp_group.add_argument('layout', metavar='Assignment sheet', type = str, widget = 'FileChooser', help = '.xlsx file with sample and assay names')
    
    args = parser.parse_args(sys.argv[1:])
    return args