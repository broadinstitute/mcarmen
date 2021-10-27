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
    
    # Hit calling
    hit_group = parser.add_argument_group("Advanced")
    hit_group.add_argument('-threshold', '--threshold', metavar = 'Threshold for hit calling', type=float, default=defaults['threshold'], help='usually between 1.2 and 2')
    hit_group.add_argument("-toi", "--toi", metavar= 'Timepoint for hit calling', type = str, default=defaults['timepoint'], help="usually after 30 to 60 min, e.g. t9")
    # controls
    hit_group.add_argument('-ectrl','--ectrl', metavar = 'Patient negative sample', type=str, default=defaults['Ectrl'])
    hit_group.add_argument('-ntcctrl', '--ntcctrl', metavar= 'Negative RT-PCR control', type = str, default=defaults['NTCctrl'])
    hit_group.add_argument('-cpcctrl','--cpcctrl', metavar = 'Pooled positive RT-PCR control',type=str,default=defaults['CPCctrl'])
    hit_group.add_argument('-ndcctrl','--ndcctrl', metavar = 'Sample MM without Mg',type=str,default=defaults['NDCctrl'])
    hit_group.add_argument('-dmctrl','--dmctrl', metavar = 'Assay MM without crRNA',type=str,default=defaults['DMctrl'])
    hit_group.add_argument('-wctrl','--wctrl', metavar = 'Water only',type=str,default=defaults['Wctrl'])
    
    #Order of Sorting
    #map_group = parser.add_argument_group('Optional')
    hit_group.add_argument("-samplelistsort", "--slsort", metavar = 'Order of samples', choices=['original', 'alphabetical'], default = defaults['slsort'])
    hit_group.add_argument("-assaylistsort", "--alsort", metavar = 'Order of assays', choices=['original', 'alphabetical'], default = defaults['alsort'])
    
    # Version control and debugging
    # map_group.add_argument('-V', '--version', action='version', version='%(prog)s {0}'.format(version)) # for command line tool if we want to offer version control
    hit_group.add_argument('-v', '--verbose', dest='verbosity', choices=['error', 'warning', 'info', 'debug'], default = defaults['verbosity'], \
                           help='Set logging level (default {verbosity})'.format(**defaults))
    
    args = parser.parse_args(sys.argv[1:])
    return args