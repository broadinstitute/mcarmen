'''

CARMEN RVP analysis v6.
For 192.192 IFCs only.

'''

# General modules
import sys
import logging
import time
from os import path, makedirs
from gooey import Gooey, GooeyParser

# import submodules
from parser import ArgParse
from analysis import ReadData, ItemsToList, LabelInputs
from plot import PlotHeatmap
from hitcalling import HitCallingFoldThreshold, CheckOutliers, CheckEC, CheckNDC, CheckDM, CheckEx, ConsiderControls, getCtrlMean


@Gooey(program_name='CARMEN RVP analysis',
      program_description = 'Analyze CARMEN RVP data from Fluidigm Biomark HD instrument.',
      tabbed_groups=True,
      language = 'english',
      image_dir = '../img/',
      clear_before_run = True,
      menu=[{
        'name': 'About',
        'items': [{
                'type': 'AboutDialog',
                'menuTitle': 'About',
                'name': 'CARMEN RVP analysis Demo',
                'description': 'An example description',
                'version': '1.0.6',
                'copyright': '2021',
                'website': 'Github Link',
                'license': 'MIT'
            }, {
                'type': 'Link',
                'menuTitle': 'Visit Our Webite',
                'url': 'https://www.sabetilab.org/'
            }]
          }]
      )

# these lines are for packaging with PyInstaller
class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)
sys.stdout = Unbuffered(sys.stdout)

# Get the version:
version = {}
with open(path.join(path.abspath(path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

# A function to quit with an error:
def error(msg, exit_code=1):
    logging.error(msg)
    sys.exit(exit_code)

def main():
    
    # Define the user defaults:
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    defaults = {'xname': timestamp, 'outdir': '', 'rawdata': '', 'layout': ''}
    
    args = ArgParse(defaults)
    args.verbosity = 'info'
    
    # Check that the output directory exists, and make it if not:
    output_dir = path.expanduser(path.abspath(args.outdir))
    if not path.exists(output_dir):
        try:
            makedirs(output_dir)
        except: 
            logging.basicConfig(level = logging.DEBUG)
            error('failed to create output directory {}'.format(output_dir))
    if not path.isdir(output_dir): 
        logging.basicConfig(level = logging.DEBUG)
        error('specified output {} is not a directory'.format(output_dir))   
    # Build the output file prefix:
    output_prefix = path.join(args.outdir, args.xname)
    
    # Set up logging based on the verbosity level set by the command line arguments:
    logfilepath = '{}logfile.log'.format(output_prefix)
    logging.basicConfig(format='%(asctime)s, %(name)s %(levelname)s %(message)s',
                        datefmt='%m-%d %H:%M',
                        level = logging.DEBUG,
                        filename= logfilepath,
                        filemode= 'w')
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(args.verbosity.upper())
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    
    logging.debug('Logfile is saved to {}'.format(logfilepath))
    logging.debug('output file prefix is {}'.format(output_prefix))
    
    # Check if everything is there
    for inputfile in [args.layout, args.rawdata]:      
        if not path.exists(inputfile):
            error('specified file {} does not exist'.format(inputfile))      
    
    # Setup of arguments that are NOT given by the user, since they are fixed
    args.toi = 't12'
    args.threshold = 1.8
    args.alsort = 'orginial'
    args.slsort = 'original'
    args.ntcctrl = 'NTC'
    args.cpcctrl = 'CPC'
    args.ectrl = 'EC'
    args.ndcctrl = 'NDC'
    args.dmctrl = 'no-crRNA'
    args.wctrl = 'water'
    
##TODO## Check if all controls are present in the layout sheet
    
    # Read in data
    logging.info("Reading in rawdata from {}".format(args.rawdata))
    try:
        probe_df, ref_df, bkgdprobe_df, bkgdref_df = ReadData(args.rawdata)
    except:
        error("Reading the raw data from {} failed.".format(args.rawdata))
    
    # Background subtraction
    probe_df = probe_df - bkgdprobe_df
    ref_df = ref_df - bkgdref_df
    # Normalization with ROX
    signal_df = probe_df/ref_df
    
    # Assign sample and assay names
    try:
        LabelInputs(signal_df, args.layout)
        logging.info("Assays and Samples have been labeled.")
    except:
        error("Unable to label assays and samples. Please check the layout sheet.")
    
    # Save signal dataframe to csv
    signal_df.to_csv("{}_signal.csv".format(output_prefix))  
    logging.info("Normalized data saved to {}_signal.csv".format(output_prefix))

    # Get median data for timepoint of interest
    x_med = signal_df.groupby(['assay', 'sample']).median()[args.toi].unstack()
    x_med.index.names=['']
    x_med.columns.names=['']
    median_df = x_med
    
    # Create a list of all assays and all samples
    a_list = ItemsToList(args.layout, 'assays', sort = args.alsort)
    s_list = ItemsToList(args.layout, 'samples', sort = args.slsort)
    
    ####################### Hit calling
    
    neg_ctrl_mean, neg_ctrl_std = getCtrlMean(median_df,args.ntcctrl,args)
    pos_ctrl_mean, pos_ctrl_std = getCtrlMean(median_df,args.cpcctrl,args)
    
    # Overall dynamic range of the assay --> turns whole assay invalid
    """
    The NTC and CPC are the ground truth for the hit calling. Since the validation of CPC and NTC depends on their group mean and standard deviation respectively, they have to be distinctively different from each other. For the assay to be valid, the ratio of separation band to dynamic range has to be greater than 0.2. The dynamic range is defined as the difference of the NTC mean and CPC mean. The separation band is the range between the CPC mean subtracted by three standard deviations of the CPCs and the NTC mean plus three standard deviations of the NTCs. 
    """
    separation_band = (pos_ctrl_mean - 3*pos_ctrl_std) - (neg_ctrl_mean + 3* neg_ctrl_std)
    dynamic_range = abs(pos_ctrl_mean - neg_ctrl_mean)
    zfactor = separation_band/dynamic_range
    
    if zfactor > 0.2:
        logging.info('Z-Factor is {}'.format(zfactor))
    else:
        logging.error("Run invalid: the dynamic range of this assay is too low (Z-factor {}).".format(zfactor))
        logging.debug("CPC mean = {} stdev = {}; NTC mean = {} stdev = {}".format(pos_ctrl_mean, pos_ctrl_std,neg_ctrl_mean, neg_ctrl_std))
        PlotHeatmap(median_df, None, args.toi, a_list, s_list, output_prefix+'_LowDynamicRange')
        logging.warning("A heatmap has been generated nevertheless, but without hit calling.")
        sys.exit()
        
    # Check for RT-PCR control outlier --> turns assay invalid
    NTC_outlier = CheckOutliers(median_df,args.ntcctrl,'negative')  
    CPC_outlier, pass_cpc = CheckOutliers(median_df,args.cpcctrl,'positive')
    
    # perform hit calling on remaining samples
    hit_df, quanthit_df = HitCallingFoldThreshold(median_df,args.threshold,a_list,s_list,args.ntcctrl)
    quanthit_df.T.to_csv("{}_HitQuantification.csv".format(output_prefix)) 
    PlotHeatmap(median_df, hit_df, args.toi, a_list, s_list, output_prefix+'_initial_raw_hitcalling')
    
    # Extraction control should be negative for all targets and positive for RNAseP
    EC_outlier, pass_ec = CheckEC(hit_df,a_list,args.ectrl)
    
    if pass_ec or pass_cpc == False:
        hit_df = hit_df.applymap(lambda x: 'invalid')
        hit_df.T.to_csv("{}_hits.csv".format(output_prefix))
        PlotHeatmap(median_df, hit_df, args.toi, a_list, s_list, output_prefix)
        logging.info("Heatmap with hits is saved as {}_heatmap_{}.png".format(output_prefix, args.toi))
        sys.exit()
        
    # Negative Detection Control (NDC) for sample mastermix should be negative for all targets in W-det-NoMg
    NDC_outlier = CheckNDC(hit_df,a_list,args.ndcctrl)
    
    # Detection control for assay mastermix should be negative for all sample in no-crRNA assay
    DM_soutlier = CheckDM(hit_df,s_list,args.dmctrl) # this outlier list is with samples instead of assays
    
    # All samples should be positive for RNAseP or at least another sample
    Ex_soutlier = CheckEx(hit_df,s_list,a_list,args) # check if extraction for the sample was successfull
    
    # annotate samples that didn't pass the control
    hit_df = ConsiderControls(hit_df,a_list,s_list,args,NTC_outlier, CPC_outlier, EC_outlier, NDC_outlier, DM_soutlier, Ex_soutlier)
    
    hit_df.T.to_csv("{}_hits.csv".format(output_prefix))            
    logging.warning("Hit calling data frame is generated and saved to {}_hits.csv.".format(output_prefix))
     
    # Plot heatmap with overlayed hits
    PlotHeatmap(median_df, hit_df, args.toi, a_list, s_list, output_prefix)
    logging.info("Heatmap with hits is generated and saved as {}_heatmap_{}.".format(output_prefix, args.toi))
    
    

############### EXECUTE #################
    
if __name__ == '__main__':
    main()

