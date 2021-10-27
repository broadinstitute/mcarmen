import pandas as pd
import numpy as np
import logging
import itertools

def HitCallingFoldThreshold(frame, fold_threshold, assaylist, samplelist, ctrl):
    # returns dataframe with hit yes or no
    data = frame[samplelist].reindex(assaylist)
    pdf = pd.DataFrame(columns = samplelist, index = assaylist)
    pdf_quant = pd.DataFrame(columns = samplelist, index = assaylist)
    for target in samplelist:
        for guide in assaylist:
            # Hitcalling
            fold = np.median(data.loc[guide, target])/np.median(frame.loc[guide, ctrl])
            pdf_quant.loc[guide,target] = fold
            if fold > fold_threshold:
                pdf.loc[guide, target] = 'positive'
            else:
                pdf.loc[guide, target] = 'negative'
    return pdf, pdf_quant


def getCtrlMean(df,ctrl,args):
    ctrl_values = df[ctrl].tolist()
    assaylist = df.index.values.tolist()
    
    if ctrl == args.cpcctrl:
        # Exclude no_crRNA control for CPC
        noCRNAidx = assaylist.index('no-crRNA')
        del ctrl_values[noCRNAidx]
    
    mean = np.mean(ctrl_values)
    std = np.std(ctrl_values)
    
    return(mean, std)

##### Process controls

def CheckOutliers(df,ctrl,direction, pass_ntc, pass_cpc):
    """ Compare negative/positive controls with each other using the Z-score.
    
The NTC consists of using molecular-grade nuclease-free water in RT-PCR reactions instead of RNA. The NTC reactions for all primer and crRNA sets should not generate any signal. An assay is defined as positive, if it falls outside of three standard deviations of the mean of all NTCs. The NTC with the no-crRNA assay is only an outlier, if the signal is more than three standard deviations above the mean. If it is lower, no error occurs, since the background signal is generally lower for the no-crRNA assay. If any of the NTC reactions (RT-PCR) generates an outlying signal, sample contamination may have occurred. This will invalidate the assay of the positive viral marker.

The Combined RVP Positive Control consists of nine components (partial synthetic nucleic acids), combined as directed in the Reagent and Controls Preparation section to make the final control sample for the reaction. Combined RVP Positive Control consists of partial synthetic targets at 1e3 copy/ul. The Combined RVP Positive Control is prepared alongside each batch of clinical specimens and should be positive for all nine targets in the CARMEN RVP Assay.
If the CARMEN RVP PC generates a negative result for any target, this indicates a possible problem with primer mixes used in RT-PCR reactions or crRNA-cas13 detection and will invalidate the assay with the negative CPC. The assay is defined as negative, if it falls outside of three standard deviations of the mean of all CPCs. The CPC with the no-crRNA assay is expected to be negative, is excluded from the CPC mean, and will be invalid if the signal falls within three standard deviations of the mean of the CPCs.
    """
    ctrl_values = df[ctrl].tolist()
    logging.info("ctrl values: {} - {}".format(ctrl_values, len(ctrl_values)))
    assaylist = df.index.values.tolist()
    logging.info("assaylist: {} - {}".format(assaylist, len(assaylist)))
    outliers = []

    threshold = 3
    mean = np.mean(ctrl_values)
    std = np.std(ctrl_values)

    for y in range(len(ctrl_values)):
        # Formula for Z score = (Observation — Mean)/Standard Deviation
        z_score= (ctrl_values[y] - mean)/std
        logging.info("ctrl: {} - assay {}".format(ctrl_values[y], assaylist[y]))
        
        # Check no crRNA and invalidate run if no-crRNA is positive for the CPC or the NTC
        if assaylist[y] == 'no-crRNA' and direction == 'positive':
            if np.abs(z_score) < threshold:
                # this is a problem, because it should be negative. The whole assay will be invalid.
                logging.warning("Run is invalid, because the CPC is not below threshold for no-crRNA (not negative enaugh).")
                pass_cpc = False    
        elif assaylist[y] == 'no-crRNA' and direction == 'negative':
            # value for no-crRNA in NTC has to be around the mean or lower than other negative samples, else:
            if ctrl_values[y] > mean + threshold * std:
                logging.warning("Run is invalid, because the NTC is above the threshold for no-crRNA")
                pass_ntc = False
        # Check that RNaseP is positive for CPC
        elif assaylist[y] == 'RNaseP' and direction == 'positive':
            if np.abs(z_score) > threshold:
                # this is a problem, because it should be positive and have a small z_score (thus similiar to the others). The whole assay will be invalid.
                logging.warning("Run is invalid, because the CPC is above the threshold for RNAseP.")
                pass_cpc = False 
        # Check that RNaseP is negative for NTC
        elif assaylist[y] == 'RNaseP' and direction == 'negative':
            # value for RNaseP in NTC has to be around the mean of all the other negative samples, else:
            if np.abs(z_score) > threshold:
                logging.warning("Run is invalid, because the NTC is positive for RNaseP")
                pass_ntc = False
        else:
            if np.abs(z_score) > threshold:
                outliers.append(assaylist[y])
                logging.debug("The {} ctrl {} with assay {} falls outside of 3 st.dev from the mean of all ctrls.".format(direction, ctrl,assaylist[y]))
            
    if direction == 'positive' and len(outliers)>0:
        logging.warning('CPC outlier {}'.format(outliers))
    elif direction == 'negative'and len(outliers)>0:
        logging.warning('NTC outlier {}'.format(outliers))
    
    return outliers, pass_ntc, pass_cpc


def CheckEC(df,assaylist,ctrl):
    """
    Extraction Negative Control (EC) is used as an RNA extraction procedural control to demonstrate successful recovery of nucleic acid, extraction reagent integrity and as a control for cross-contamination during the extraction process. 
    
    The EC control consists of a confirmed negative patient sample. Purified nucleic acid from the EC should yield a positive result with the RP primer and crRNA set and negative results for viral specific targets. 
    
    If the EC generates a negative result for RP, this indicates a potential problem with the extraction process. This will invalidate the whole run. Thus, repeat extraction, RT-PCR and crRNA-Cas13 testing for all specimens that had been extracted alongside the failed EC. If EC generates positive results for any of the viral markers, this may be an indication of possible cross-contamination during extraction, RT-PCR or crRNA-Cas13 reaction set-up. In this case, only the assay of the positive viral marker will be turned invalid.
    """
    
    outliers = []
    passed = True
    
    for guide in assaylist:
        if guide == 'RNaseP':
            if df.loc[guide, ctrl] != 'positive':
                logging.debug('guide: {}, value {}'.format(guide,df.loc[guide, ctrl]))
                logging.warning('EC is not positive for RNaseP and run is invalid because of failed extraction control.')
                passed = False
        elif guide == 'no-crRNA':
            if df.loc[guide, ctrl] == 'positive':
                logging.debug('guide: {}, value {}'.format(guide,df.loc[guide, ctrl]))
                logging.warning('EC is positive for no-crRNA and run is invalid.')
                passed = False
        else:
            if df.loc[guide, ctrl] != 'negative':
                logging.debug('EC is not negative for {}'.format(guide))
                outliers.append(guide)
    
    logging.warning('EC not successfull for {}'.format(outliers)) 
    return outliers, passed
    
def CheckNDC(df,assaylist,ctrl):
    """
    The negative detection controls (NDC) consist of using molecular-grade nuclease-free water in crRNA-Cas13 detection reactions instead of RNA without the presence of Magnesium (Mg++)in the reaction mastermix. If the negative detection control is positive for a viral target, all samples of this assay will be invalid.
    """
    outliers = []
    pass_ndc = True
    
    for guide in assaylist:
        if guide == 'RNaseP':
            if df.loc[guide, ctrl] == 'positive':
                logging.debug('guide: {}, value {}'.format(guide,df.loc[guide, ctrl]))
                logging.warning('NDC is positive for RNaseP and run is invalid.')
                pass_ndc = False
        elif guide == 'no-crRNA':
            if df.loc[guide, ctrl] == 'positive':
                logging.debug('guide: {}, value {}'.format(guide,df.loc[guide, ctrl]))
                logging.warning('NDC is positive for no-crRNA and run is invalid.')
                pass_ndc = False            
        elif df.loc[guide,ctrl] != 'negative':
            logging.debug('Sample mastermix without Mg is not negative for {} in {}'.format(guide,ctrl)) 
            outliers.append(guide)
            
    logging.warning('NDC not successfull for {}'.format(outliers))        
    return outliers, pass_ndc


##### Internal controls


def CheckDM(df,samplelist,ctrl):
    """
    The second negative control (no crRNA) consists of a detection master mix with molecular-grade nuclease-free water instead of crRNA. If a sample is not negative no the no crRNA control, all assays for this sample will be invalid. 
    """
    outliers = []
    for sample in samplelist:
        if df.loc[ctrl,sample] != 'negative':
            logging.debug('Detection Assay MM control is not negative in sample {} for {} assay'.format(sample,ctrl))
            outliers.append(sample)       
    logging.warning('Detection Assay MM control (no crRNA) is not negative for {}'.format(outliers))
    return outliers
    
def CheckEx(df,samplelist,assaylist,args):
    """ Sample-specific extraction control. 
    For all samples, at least one assay has to be positive otherwise this indicates a problem with extraction and the sample will be invalid. Repeat extraction, RT-PCR and crRNA-Cas13 testing for this specimen.
    """
    outliers = []
    dontcheck = [args.ectrl,args.ntcctrl,args.cpcctrl,args.ndcctrl,'water_exclude',args.wctrl] # exclude controls from this analysis
    for sample in samplelist:
        if sample in dontcheck:
            continue
        count1 = sum(1 for x in df[sample] if x == 'positive')
        if count1 == 0:
            outliers.append(sample)
            
    logging.warning('Sample specific extraction control (RNaseP or positive for at least one other viral target) not successfull for {}'.format(outliers))
    return outliers
        
##### Finalize hit calling

def ConsiderControls(df,assaylist,samplelist,args,NTCout,CPCout,ECout,DSout,DMout,Exout):
    """ Replace called hit (1) or no hit (0) with invalid result (-1), or causing invalid result (-2) in hit dataframe.
    """
    resultx = 'invalid' # invalid
    resulty = 'invalid' # the sample causing the invalid result. For clinical reporting, both are named the same
    
    #allguidesoutliers = NTCout + CPCout + ECout + DSout
    allguidesoutliers = list(itertools.chain(NTCout, CPCout,ECout,DSout))
    logging.info('Invalid guides: {}'.format(allguidesoutliers[0]))
    allsamplesoutliers = DMout + Exout
    logging.info('Invalid samples: {}'.format(allsamplesoutliers))
    for guide in assaylist:
        for sample in samplelist:
            # skip if guide/sample combination is fine
            if guide not in allguidesoutliers and sample not in allsamplesoutliers:
                continue
            # If negative RT-PCR control invalid, all other samples are invalid
            if guide in NTCout:
                if sample == args.ntcctrl:
                    df.loc[guide, sample] = resulty # this is causing the invalid result
                elif sample == args.cpcctrl and guide in CPCout:
                    df.loc[guide, sample] = resulty # this is also causing the invalid result
                elif sample == args.cpcctrl:
                    continue
                else:
                    df.loc[guide, sample] = resultx # all the others are invalid
            # If positive RT-PCR control invalid, all other samples, besides the negative control are invalid
            elif guide in CPCout:
                if sample == args.ntcctrl:
                    continue
                elif sample == args.cpcctrl:
                    df.loc[guide, sample] = resulty # this is causing the invalid result
                else:
                    df.loc[guide, sample] = resultx # all the others are invalid      
            # if any other control is invalid, all the samples besides controls are invalid
            else:
                if guide in ECout:
                    if sample == args.ectrl:
                        df.loc[guide, sample] = resulty # this is causing the invalid result
                    elif sample in [args.ntcctrl,args.cpcctrl,args.ndcctrl,args.dmctrl]:
                        continue # this sample represents the ground truth and shouldn't be changed
                    else:
                        df.loc[guide, sample] = resultx # invalid
                if guide in DSout:
                    if sample == args.ndcctrl:
                        df.loc[guide, sample] = resulty # this is causing the invalid result
                    elif sample in [args.ntcctrl,args.cpcctrl,args.ectrl,args.ndcctrl,args.dmctrl]:
                        continue # this sample represents the ground truth and shouldn't be changed
                    else:
                        df.loc[guide, sample] = resultx # invalid
                if sample in DMout:
                    if guide == args.dmctrl:
                        df.loc[guide, sample] = resulty # this is causing the invalid result
                    else:
                        df.loc[guide, sample] = resultx # invalid
                if sample in Exout:
                    if guide == 'RNAseP' or 'RNaseP':
                        df.loc[guide, sample] = resulty # this is causing the invalid result
                    elif guide == args.dmctrl:
                        continue
                    else:
                        df.loc[guide, sample] = resultx # invalid
    return df

                    
                    