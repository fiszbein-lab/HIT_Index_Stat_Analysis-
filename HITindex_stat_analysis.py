
import sys
import pandas as pd
import numpy as np
import collections
import math
import scipy
import statsmodels.api as sm
import argparse
import statsmodels.stats as sms


def process_nup_ndown(file):
    """ read in the HITindex data and return a dataframe that contain difference
    between nup and ndown
    """
    if file[-6:] == 'AFEPSI':
        position = "start"
    else:
        position = "end"
    # determine AFE or ALE

    data = pd.read_csv(file, sep="\t")

    #compute the difference
    if position == "start":
        data["nup_ndown_difference"] = data['nDOWN'] - data['nUP']
    else:
        data["nup_ndown_difference"] = data['nUP'] - data['nDOWN']
    return data


def dictionary_build(df, dic):
    """input the value into a dictionary for further process, this function will
    be repeatly called until all samples are entered into the dic.
    dic is a nested dic created with collections.defaultdict(lambda:collections.defaultdict(dict))
    """

    for index, row in df.iterrows():
        current_gene = row['gene']
        current_exon = row['exon']
        if current_gene not in dic:
            dic[current_gene][current_exon]['read'] = [int(row['nup_ndown_difference'])]
        elif current_exon not in dic[current_gene]:
            dic[current_gene][current_exon]['read'] = [int(row['nup_ndown_difference'])]
        else:
            dic[current_gene][current_exon]['read'] += [int(row['nup_ndown_difference'])]
    return dic




def merge_samples(proceed_samples,gene_dic, samples_name, outname = 'merged_samples'):
    """this function merge all the samples into a single file
    result will be a dataframe file contain all those informations
    proceed_samples should be a list contain dataframe of samples
    samples_name is a list that contain the name of the samples
    """
    sample_dictionaries = []
    output = open(outname, 'w')
    Type = samples_name[1][-6:]


    first_line = 'gene exon'
    for name in samples_name:
        first_line += f' {name[:-7]}.difference'
    for name in samples_name:
        first_line += f' {name[:-7]}.PSI'
    first_line += ' type'

    print(first_line, file = output)


    for sample in proceed_samples:
        # sample is a dataframe
        # read in all value of samples into dics, though it is not necessary
        tem_dic = collections.defaultdict(dict)
        for index, row in sample.iterrows():
            current_gene = row['gene']
            current_exon = row['exon']
            current_value = row['nup_ndown_difference']
            current__PSI_value = row[f'{Type}']
            tem_dic[current_gene][current_exon] = [current_value, current__PSI_value]

        sample_dictionaries += [tem_dic]

    for gene in gene_dic:
        for exon in gene_dic[gene]:
            current_str = f'{gene} {exon}'
            for dic in sample_dictionaries:
                if gene not in dic or exon not in dic[gene]:
                    current_str += ' 0'
                else:
                    current_str += f' {dic[gene][exon][0]}'

            for dic in sample_dictionaries:
                if gene not in dic or exon not in dic[gene]:
                    current_str += ' 0'
                else:
                    current_str += f' {dic[gene][exon][1]}'
            if samples_name[0][-6:] == "ALEPSI":
                current_str += ' ALE'
            else:
                current_str += ' AFE'



            print(current_str, file = output)

    # must close the output first or the last several lines of data will in the
    #cache and causing the lack of several lines in the new dataframe
    output.close()
    result = pd.read_csv(outname, sep=" ")

    # return a dataframe that can be use for further process
    return(result)

def merge_AFEPSI_ALEPSI(df_AFE,df_ALE,outname = 'merged_samples'):
    """this function merge AFEPSI and ALEPSI for further computation
    """
    result = pd.concat([df_AFE,df_ALE])
    result.to_csv(outname, sep= ' ', index=False)
    return result


def built_gene_exon_dic(sample_df):
    """ this function built a dictionary that contain all gene and exons in samples
    """
    result_dic = collections.defaultdict(dict)
    for index, row in sample_df.iterrows():
        current_gene = row['gene']
        current_exon = row['exon']
        result_dic[current_gene][current_exon] = str(row['type'])

    return result_dic



def calculate_likehood(df, list_of_exon, list_of_sample_name, current_exon, threshold, exon_threshold):
    """this function is a helper function of model_computation, it calculate
    the likehood of models
    """
    df_test = df.loc[df['condition'] == 1]
    df_control = df.loc[df['condition'] == 0]
    test_count = df_test.shape[0]
    number_of_values = df.shape[0]



    current_exon_lines = df.loc[df[current_exon] == 1]
    current_exon_count = current_exon_lines['value'].sum()

    if df['value'].sum() < threshold or df_test['value'].sum() == 0\
        or df_control['value'].sum() == 0 or current_exon_count < exon_threshold or test_count < 2:
        # mark untestable condition
        p_val = -1

    else:

        # full model
        y1 = df['log_value']

        x1_varaible =  list_of_exon + list_of_sample_name + ['condition_exon']
        x1 = df[x1_varaible]
        x1 = sm.add_constant(x1)
        full_model = sm.OLS(y1, x1).fit()
        # calculate the log-likehood test
        full_ll = full_model.llf



        #reduced model
        y2 = df['log_value']
        x2_varaible =  list_of_exon + list_of_sample_name
        x2 = df[x2_varaible]
        x2 = sm.add_constant(x2)
        reduced_model = sm.OLS(y2, x2).fit()
        reduced_ll = reduced_model.llf
        LR_statistic = -2*(reduced_ll - full_ll)
        p_val = scipy.stats.chi2.sf(LR_statistic,1)

    #print(p_val)
    return p_val


def PSI_value_dic_construct(samples_name, merged_samples_df):
    """this function return a dic that contain values for PSI values for samples
    """
    PSI_dic = collections.defaultdict(lambda:collections.defaultdict(dict))
    for index, row in merged_samples_df.iterrows():
        current_gene = row['gene']
        current_exon = row['exon']
        for i in samples_name:
            current_PSI_value = row[f'{i}.PSI']
            PSI_dic[current_gene][current_exon][i] = current_PSI_value

    return PSI_dic

def outlier_detection(sample_type, df_row, outlier_treatment = 'default', outlier_method = 'cooks'):
    """this is the function use to detect outliar, it will take in a list of samples
    names and a df_row that contain the read information about a specific exon and check
    whether there is an outliar. Due to the sample size limitation, test and control group
    will be merge together for the computation if each group have fewer than 5 samples.
    The function will output an list of sample_names with no outliar for next step compuation
    and a list of outliar. It use 1.5 IQR to to see whether there is an outliar
    returned result in format: [usable_sample_list, outlier]
    """

    # to check whether both condition have at least five samples each(not in use right now)
    test_sample = 0
    control_sample = 0
    test_sample_list = []
    control_sample_list = []
    for i in sample_type:
        if i[1] == "control":
            control_sample += 1
            control_sample_list += [i[0]]
        else:
            test_sample += 1
            test_sample_list += [i[0]]
    list_of_sample_name = test_sample_list + control_sample_list

    if outlier_method == 'cooks':
        result = cooks_outlier_calculation(list_of_sample_name, sample_type, df_row, 'merged', outlier_treatment)
    elif outlier_method == 'none':
        result = [list_of_sample_name, []]
    else:
        if outlier_treatment == "merge":
            result = iqr_outlier_calculation(list_of_sample_name, sample_type, df_row, 'merged')
        elif outlier_treatment == 'separate':
            test_result = iqr_outlier_calculation(test_sample_list, sample_type, df_row, 'separated')
            control_result = iqr_outlier_calculation(control_sample_list,  sample_type, df_row, 'separated')
            result = [test_result[0]+control_result[0],test_result[1]+control_result[1] ]
        elif (len(test_sample_list) < 5 and len(control_sample_list) < 5):
            result = iqr_outlier_calculation(list_of_sample_name,  sample_type, df_row, 'merged')
        else:
            test_result = iqr_outlier_calculation(test_sample_list,df_row, 'separated')
            control_result = iqr_outlier_calculation(control_sample_list,  sample_type, df_row, 'separated')
            result = [test_result[0]+control_result[0],test_result[1]+control_result[1] ]

    return result


def iqr_outlier_calculation(sample_names,  sample_type, df_row, data_type):
    """this is the helper function for outlier_dectection, which do all the conputation
    this function will return a list of list that contain information about usable samples
    and outliers
    uses iqr method
    """

    list_of_value = []

    for i in sample_names:
        list_of_value += [float(df_row[f'{i}.PSI'])]

    array_of_values = np.array(list_of_value)
    if data_type == 'merged':
        q1 =  np.percentile(array_of_values, 25, interpolation = 'lower')
        q3 =  np.percentile(array_of_values, 75, interpolation = 'higher')
    else:
        q1 =  np.percentile(array_of_values, 25)
        q3 =  np.percentile(array_of_values, 75)
    IQR = q3 - q1
    outlier_low_limit = q1 - (1.5*IQR)
    outlier_high_limit = q3 + (1.5*IQR)


    usable_samples = []
    outlier = []
    for i in sample_names:
        current_PSI_value = float(df_row[f'{i}.PSI'])
        if current_PSI_value < outlier_low_limit or current_PSI_value > outlier_high_limit:
            outlier += [i]
        else:
            usable_samples += [i]
    return [usable_samples, outlier]

def cooks_outlier_calculation(sample_names, sample_type, df_row, data_type, outlier_treatment):
    """this is a helper function for outlier_detection, which doees the computation
    this funciton will return a list of lists that contain information about usable samples
    and outlier
    uses cooks method
    """
    list_of_value = []
    samp = []
    for i in sample_type:
        list_of_value += [float(df_row[f'{i[0]}.PSI'])]
        if i[1] == "control":
            samp += [0]
        else:
            samp += [1]
    Y = list_of_value
    X = samp
    X = sm.add_constant(X)
    model = sm.OLS(Y, X).fit()
    np.set_printoptions(suppress=True)
    influence = model.get_influence()
    cooks = influence.cooks_distance
    if outlier_treatment == "4/n":
        usable_samples = [sample_names[i] for i in [i for i,v in enumerate(cooks[0]) if v < 4.0/len(sample_names)]]
        outlier = [sample_names[i] for i in [i for i,v in enumerate(cooks[0]) if v >= 4.0/len(sample_names)]]
    elif outlier_treatment == "4*mean":
        usable_samples = [sample_names[i] for i in [i for i,v in enumerate(cooks[0]) if v < 4.0*np.mean(cooks[0])]]
        outlier = [sample_names[i] for i in [i for i,v in enumerate(cooks[0]) if v >= 4.0*np.mean(cooks[0])]]
    elif outlier_treatment == "1":
        usable_samples = [sample_names[i] for i in [i for i,v in enumerate(cooks[0]) if v < 1]]
        outlier = [sample_names[i] for i in [i for i,v in enumerate(cooks[0]) if v >= 1]]
    else:
        print("ERROR: must select proper outlierDetection for the selected outlierMethod (4/n or 4*mean or 1 for cooks, default or none or merge or separate for iqr)")
    return [usable_samples, outlier]

def less_than_half_expression_detection(df, usable_samples, sample_type):
    """this function check whether less than half of the control or test samples
    express this gene
    """
    num_control = 0
    num_test = 0
    num_zero_control = 0
    num_zero_test = 0
    dic_sample = {'test': [], 'control': []}
    for i in sample_type:
            dic_sample[i[1]] += [i[0]]
    for sample in usable_samples:
        if sample in dic_sample['test']:
            num_test += 1
            temp_df = df.loc[df[sample] == 1]
            if temp_df['value'].sum() == 0:
                num_zero_test += 1
        else:
            num_control += 1
            temp_df = df.loc[df[sample] == 1]
            if temp_df['value'].sum() == 0:
                num_zero_control += 1
    if num_test < 2*num_zero_test or num_control < 2*num_zero_control:
        return True
    else:
        return False



def model_computation(gene_exon_dic, merged_result, sample_type, PSI_dic,\
                      output_name = 'Hitindex_stat_interpretation', threshold = 0, \
                          exon_threshold = 0, critical_value = '0.01', outlier_treatment = 'default', outlier_method = 'cooks', multiple_testing = 'fdr_bh'):
    """this is the function that do the computation. It will built linear regression
    model based on exon in each gene. It require input of gene_exon_dic that is a dic
    that contain all gene and exon. merged_result will be result from previous steps
    sample_type will be a list that contain lists of information about conditions: for example:
    [['sample_1', 'test'], ['sample_2', 'control']]
    """

    result_dic = collections.defaultdict(dict)
    outlier_dic = collections.defaultdict(dict)
    p_vals = []
    # the structure of this dic will be gene:exon:[p_value, a=0.1, a=0.05, a=0.01]
    # boolean values for different threshold
    for gene in gene_exon_dic:
        df = merged_result[merged_result['gene'] == gene]
        for exon in gene_exon_dic[gene]:
            current_exon_type = gene_exon_dic[gene][exon]
            #print(f'running for {current_exon_type} {exon} in {gene}')
            # computate two models and perfrom likehood test for all exon
            list_of_exon = []
            dataframe_string = {'value':[],'log_value': [], 'condition':[], 'condition_exon': [], 'ALE': [], 'AFE': []}
            for exon_name in gene_exon_dic[gene]:
                list_of_exon += [exon_name]
                dataframe_string[exon_name] = []
            outlier_result = []
            for index,row in df.iterrows():
                if row['exon'] == exon:
                    outlier_result = outlier_detection(sample_type, row, outlier_treatment, outlier_method)
                    usable_sample_names = outlier_result[0]
                    outlier_dic[row['gene']][row['exon']] = outlier_result[1]
                    #print(outlier_result[1])
                    usable_samples = []
                    for i in sample_type:
                        if i[0] in outlier_result[0]:
                            usable_samples += [i]


            for sample in usable_sample_names:
                dataframe_string[sample] = []
            # built a dictionary that can be used for the computation

            # df here contain all reads of exons assocaite with this gene
            current_exon_read_test = 0
            total_sample_read_test = 0
            current_exon_read_control = 0
            total_sample_read_control = 0
            number_ALE = 0
            number_AFE = 0
            for index, row in df.iterrows():
                #print(row)

                total_current_read = 0


                for sample in usable_samples:

                    if row['type'] == 'ALE':
                        dataframe_string['ALE'] += [1]
                        dataframe_string['AFE'] += [0]
                    else:
                        dataframe_string['AFE'] += [1]
                        dataframe_string['ALE'] += [0]


                    current_sample = sample[0]
                    current_condition = sample[1]

                    total_current_read += float(row[f'{current_sample}.difference'])

                    dataframe_string['value'] += [row[f'{current_sample}.difference']]
                    current_log_value = math.log(float(row[f'{current_sample}.difference']+1))
                    # log(x+1) to deal with 0 count value
                    dataframe_string['log_value'] += [current_log_value]

                    #control = 0, test = 1
                    if current_condition == 'test':
                        dataframe_string['condition'] += [1]
                    else:
                        dataframe_string['condition'] += [0]

                    for exon_name in list_of_exon:
                        # compute the dummy variable for the exons
                        if row['exon'] == exon_name:
                            dataframe_string[exon_name] += [1]

                        else:
                            dataframe_string[exon_name] += [0]

                    #make a dummy variable for sample
                    for sample_name in usable_sample_names:
                        if sample_name == current_sample:
                            dataframe_string[sample_name] += [1]
                        else:
                            dataframe_string[sample_name] += [0]



                    # whether it is current exon and test sample
                    if current_condition == 'test' and row['exon'] == exon:
                        dataframe_string['condition_exon'] += [1]
                        current_exon_read_test += row[f'{current_sample}.difference']
                        total_sample_read_test += 1
                    else:
                        dataframe_string['condition_exon'] += [0]

                    if current_condition == 'control' and row['exon'] == exon:
                        current_exon_read_control += row[f'{current_sample}.difference']
                        total_sample_read_control += 1

                #count the ALE and AFE for this gene and discount the exon with only outlier have expression.
                if row['type'] == 'ALE' and total_current_read > 0 :
                    number_ALE += 1
                elif row['type'] == 'AFE' and total_current_read > 0:
                    number_AFE += 1


            total_current_exon_read =  current_exon_read_test + current_exon_read_control
            average_read_test = current_exon_read_test/total_sample_read_test
            average_read_control = current_exon_read_control/total_sample_read_control


            current_p_value = 0
            computation_df = pd.DataFrame(dataframe_string)
            if number_ALE <= 1 and current_exon_type == 'ALE':
                current_p_value = -1
            elif number_AFE <= 1 and current_exon_type == 'AFE':
                current_p_value = -1
            elif total_current_exon_read == 0:
                current_p_value = -1
            elif less_than_half_expression_detection(computation_df,usable_sample_names,sample_type):
                current_p_value = -1
            else:
                current_p_value = calculate_likehood(computation_df,list_of_exon, usable_sample_names, exon, threshold, exon_threshold)


            #print('one exon finished')
            if current_p_value == -1:
                a_cust = 'untestable'
                a_005 = 'untestable'
                a_001 = 'untestable'
            else:
                a_cust = (current_p_value <= float(critical_value))
                a_005 = (current_p_value <= 0.05)
                a_001 = (current_p_value <= 0.01)
            result_dic[gene][exon] =  [current_p_value, a_cust, a_005, a_001, \
                                           current_exon_read_test, current_exon_read_control,\
                                               average_read_test, average_read_control]
            p_vals.append(current_p_value)


    a01 = []
    a05 = []
    acust = []
    adj_p_vals = multiple_testing_correction(p_vals, multiple_testing)
    for x in adj_p_vals:
        if x == -1:
            acust.append('untestable')
            a05.append('untestable')
            a01.append('untestable')
        else:
            acust.append((x <= float(critical_value)))
            a05.append((x <= 0.05))
            a01.append((x <= 0.01))


    output_file = open(f'{output_name}', 'w')
    head_string = f'gene exon p_value {critical_value} 0.05 0.01 total_condition1_read total_condition2_read average_condition1_read average_condition2_read'
    # include the PSI value for the samples which will allow a easier analysis
    for i in sample_type:
        head_string += f' {i[0]}.PSI'
    head_string += ' type'
    head_string += ' outlier'
    print(head_string, file = output_file)



    j = 0
    for gene in result_dic:
        for exon in result_dic[gene]:
            values = result_dic[gene][exon]
            string = f'{gene} {exon} {adj_p_vals[j]} {acust[j]} {a05[j]} {a01[j]} {values[4]} {values[5]} {values[6]} {values[7]}'
            for i in sample_type:
                # discard the read that ALE in one condition and AFE in others
                if gene_exon_dic[gene][exon] == 'mixed':
                    string += ' not_evaluable'
                else:
                    sample = i[0]
                    string += f' {PSI_dic[gene][exon][sample]}'
            string += f' {gene_exon_dic[gene][exon]}'



            outlier = ''
            if len(outlier_dic[gene][exon])== 0:
                outlier = 'none'
            else:
                for i in outlier_dic[gene][exon]:
                    outlier += f'{i},'
                outlier = outlier[:-1]
            string += f' {outlier}'
            j += 1
            print(string, file = output_file)



    output_file.close()
    result = pd.read_csv(f'{output_name}', sep=" ")


    # return a dataframe that can be use for further process
    return(result)


def multiple_testing_correction(p_val_list, multiple_testing_method):
    """this function performs Benjamini-Hochberg correction
    for multiple testing in the p-values, producing FDR instead of p-values
    """
    adj_p = sms.multitest.multipletests(p_val_list, method=multiple_testing_method)[1]
    readj_p = [-1.0 if i < 0 else i for i in adj_p]
    return readj_p

def biological_significant(analysis_result, sample_type, value = 0.1):
    """this function is used to find the exon that is both statiscally and biologically
    significantly expressed
    """
    result_df = pd.read_csv(analysis_result, sep= " ")
    control_sample = []
    test_sample = []
    for i in sample_type:
        if i[1] == 'test':
            test_sample += [i[0]]
        else:
            control_sample += [i[0]]

    result_df['delta_PSI'] = -1.00
    result_df['log2fc'] = 0.00
    result_df['bio_significant'] = "False"


    for index, row in result_df.iterrows():
        current_control_PSI = 0
        current_test_PSI = 0
        total_control_sample = 0
        total_test_sample = 0
        # evade the problem associate with not_evaluable data
        if row[f'{control_sample[0]}.PSI'] != 'not_evaluable':
            for i in control_sample:
                if i not in str(row['outlier']):
                    current_control_PSI += float(row[f'{i}.PSI'])
                    total_control_sample += 1
            for i in test_sample:
                if i not in str(row['outlier']):
                    current_test_PSI += float(row[f'{i}.PSI'])
                    total_test_sample += 1

            average_control_PSI = current_control_PSI/total_control_sample
            average_test_PSI = current_test_PSI/total_test_sample
            if abs(average_control_PSI - average_test_PSI) >= value:
                result_df.at[index,'bio_significant'] = "True"
            result_df.at[index,'delta_PSI'] = (average_control_PSI - average_test_PSI)
            result_df.at[index, 'log2fc'] = (math.log2((average_control_PSI+0.00001)/(average_test_PSI+0.00001)))

    result_df.to_csv(analysis_result, sep = ' ', index=False)
    return result_df



def differential_gene_extraction(input_df, outname, critical_value = '0.01'):
    """this function is use to extract gene with differntial AFE or ALE usage
    """

    total_gene_diff = 0
    total_ALE_diff = 0
    total_AFE_diff = 0
    gene_dic = collections.defaultdict(dict)
    for index, row in input_df.iterrows():
        if (row[f'{critical_value}'] == 'True') and row['bio_significant'] == 'True':
                gene_dic[row['gene']][row['exon']] = row['type']


    output = open(f'{outname}.genelist', 'w')
    print('differentially_expressed_gene exons types', file = output )
    for gene in gene_dic:
        current_exons = ''
        current_types = ''
        for exon in gene_dic[gene]:
            current_exons += f'{exon},'
            temp_type = gene_dic[gene][exon]
            current_types += f'{temp_type},'
        current_type_list = current_types.split(',')

        if 'ALE' in current_type_list:
            total_ALE_diff += 1
        if 'AFE' in current_type_list:
            total_AFE_diff += 1

        print(f'{gene} {current_exons[:-1]} {current_types[:-1]}', file = output)
        total_gene_diff += 1

    output.close()
    print('There are', total_AFE_diff, 'gene expressed differentially for AFE and', total_ALE_diff, 'for ALE in', outname)
    print('In total:', total_gene_diff)



def concise_exon_infor_extraction(input_df, outname, critical_value = '0.01'):
    """this function extract basic information about exon that differentially used as AFE or ALE between conditions. All information in the
    output file also in the basic result file, this function intent to generate a file that only contain necessary information for graph making
    and data analysis.
    """
    output = open(f'{outname}.diffexons', 'w')
    print('gene exon p-value delta_PSI log2fc type', file = output)

    for index, row in input_df.iterrows():
        if (row[f'{critical_value}'] == 'True') and row['bio_significant'] == 'True':
            print(row['gene'], row['exon'], row['p_value'], row['delta_PSI'], row['log2fc'],row['type'], file = output)

    output.close()





def stat_analysis_main_function(sample_type,  outname = 'Hitindex_stat_interpretation', biosignificant_value = 0.1,
                            threshold = 10, exon_threshold = 5, critical_value = '0.01',  outlier_treatment = 'default', outlier_method = 'cooks',
                            multiple_testing = 'fdr_bh'
                            ):
    """ this is the main function for Hitindex stat interpretation
    for easy testing purpose
    samples is a list that contain the name for the samples
    sample_type is a list that contain the name and the sample type information
    should be modified after test finished
    """

    AFE_names = []
    ALE_names = []
    names = []

    for i in sample_type:
        AFE_names += [f'{i[0]}.AFEPSI']
        ALE_names += [f'{i[0]}.ALEPSI']
        names += [i[0]]
    # built a list of name for further computation


    processing_ALE_samples = []
    for i in ALE_names:
        processing_ALE_samples += [process_nup_ndown(i)]
    # calculate the difference of nup and ndown in AFE samples

    processing_AFE_samples = []
    for i in AFE_names:
        processing_AFE_samples += [process_nup_ndown(i)]
    # calculate the difference of nup and ndown in AFE samples


    dic_ALE = collections.defaultdict(lambda:collections.defaultdict(dict))
    dic_AFE = collections.defaultdict(lambda:collections.defaultdict(dict))


    for i in processing_ALE_samples:
        dictionary_build(i,dic_ALE)
    # compute the information into a dic
    for i in processing_AFE_samples:
        dictionary_build(i,dic_AFE)
    print('Dictionary built')



    merge_AFE = merge_samples(processing_AFE_samples, dic_AFE, AFE_names, f'{outname}.mergedAFE')
    merge_ALE = merge_samples(processing_ALE_samples, dic_ALE, ALE_names, f'{outname}.mergedALE')


    print('Samples merged')

    gene_exon_dic_AFE = built_gene_exon_dic(merge_AFE)
    gene_exon_dic_ALE = built_gene_exon_dic(merge_ALE)

    print('Gene exon dics built')

    PSI_dic_AFE = PSI_value_dic_construct(names, merge_AFE)
    PSI_dic_ALE = PSI_value_dic_construct(names, merge_ALE)

    print('PSI dics built')
    print('Runing for stat analysis')

    #result = model_computation(gene_exon_dic,normalize_result,sample_type)
    print(f"{multiple_testing} method selected to correct for multiple testing")
    print(f"{outlier_treatment} selected for outlier treatment")
    print(f"{outlier_method} method selected to account for outlier/extreme PSI values")
    result_AFE = model_computation(gene_exon_dic_AFE,merge_AFE,sample_type, PSI_dic_AFE, outname, threshold, exon_threshold, critical_value, outlier_treatment, outlier_method, multiple_testing)
    result_ALE = model_computation(gene_exon_dic_ALE,merge_ALE,sample_type, PSI_dic_ALE, outname, threshold, exon_threshold, critical_value, outlier_treatment, outlier_method, multiple_testing)
    print('Stat computation finished')

    print('Start gene extraction')
    merged_result = merge_AFEPSI_ALEPSI(result_AFE, result_ALE, outname)

    #biological_significantfunction could be improved, next step should just input a df
    merged_result = biological_significant(outname, sample_type, biosignificant_value)

    differential_gene_extraction(merged_result,outname, critical_value)
    print('Differentially expressed gene extracted')
    concise_exon_infor_extraction(merged_result,outname, critical_value)
    print('Brief diffrentially expressed exons information extracted')

    return merged_result




if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Statistical analysis of HITindex results', add_help = True)


    parser.add_argument('--statAnalysis', action='store_true', default=False, help='Using Linear regression to detect differential AFE and ALE usage based on HITindex result', required = False)
    # more functionality will be added next, espectially those helper function




    # stat analysis information
    Stat_analysis_infor = parser.add_argument_group('stat analysis information')
    Stat_analysis_infor.add_argument('--condition1', type = str, metavar='', help = 'base namess of HIT index results of type1, names should be connected with ,', required=False, default="None")
    # example: 4-1744_T35_FN2.1wt_day3.5_mesoderm_progenitor_HITindex&&&5-1748_T35_FN2.1wt_day3.5_mesoderm_progenitor_HITindex
    Stat_analysis_infor.add_argument('--condition2', type = str, metavar='', help = 'base namess of HIT index results of type2, names should be connected with ,', required=False, default="None")
    Stat_analysis_infor.add_argument('--biosignificant', type = float, metavar='', help = 'threshold of biosignificancy', required=False, default=0.2)
    Stat_analysis_infor.add_argument('--output', type = str, metavar='', help = 'output name', required=False, default="Hitindex_stat_interpretation")
    Stat_analysis_infor.add_argument('--minimalGeneCount', type = int, metavar='', help = 'threshold of minimal reads that needed for a gene to be tested', required=False, default=20)
    # include this term to limit significant result cause by chances when read count of a gene is small
    Stat_analysis_infor.add_argument('--minimalExonCount', type = int, metavar='', help = 'threshold of minimal reads that needed for a exon to be tested', required=False, default=15)
    Stat_analysis_infor.add_argument('--criticalValue', type = str, metavar='', help = 'threshold of statistic significant, can be a float from 0.0 to 1.0', required=False, default= '0.01')
    Stat_analysis_infor.add_argument('--multipleTesting', type = str, metavar='', help = 'ways to correct for multiple testing, can be any input for statsmodels.stats.multitest.multipletests(), recomended: fdr_bh or bonferroni', required=False, default= 'fdr_bh')
    Stat_analysis_infor.add_argument('--outlierDetection', type = str, metavar='', help = 'ways to detect outlier, can be for iqr: default, merge, separate, can be for cooks: 4/n, 4*mean, or 1', required=False, default= '4/n')
    Stat_analysis_infor.add_argument('--outlierMethod', type = str, metavar='', help = 'method selected to detect outlier, can be cooks or iqr or none', required=False, default= 'cooks')
    args = parser.parse_args()


    if not args.statAnalysis:
        sys.exit("ERROR! Need to include at least one function: --statAnalysis")

    if args.statAnalysis:
        sample_types = []
        if args.condition1 == 'None' or args.condition2 == 'None':
            sys.exit("ERROR! Please include at least two samples for each conditions")
        if args.biosignificant < 0 or args.biosignificant > 1:
            sys.exit("ERROR! Biosignificant cannot smaller than 0 and greater than 1")
        if args.minimalGeneCount < 0 or args.minimalExonCount < 0:
            sys.exit("ERROR! minimalGeneCount and minimalExonCount cannot be negative")
        if args.outlierDetection != 'default' and args.outlierDetection != 'merge' and args.outlierDetection != 'separate'and args.outlierDetection != '4/n' and args.outlierDetection != '4*mean' and args.outlierDetection != '1':
            sys.exit('ERROR! input for outlierDetection cannot be recognized, please check input again')
        if args.outlierMethod != 'cooks' and args.outlierMethod != 'iqr' and args.outlierMethod != 'none':
            sys.exit('ERROR! input for outlierMethod cannot be recognized, please check input again')
        if args.outlierMethod == 'cooks' and (args.outlierDetection == 'default' or args.outlierDetection == 'merge' or args.outlierDetection == 'separate'):
            sys.exit('ERROR! You are using cooks outlier detection method, please only use default, merge, or separate for iqr')
        if args.outlierMethod == 'iqr' and (args.outlierDetection == '4/n' or args.outlierDetection == '4*mean' or args.outlierDetection == '1'):
            sys.exit('ERROR! You are using iqr outlier detection method, please only use 4/n, 4*mean, or 1 for cooks')
        if args.multipleTesting != 'bonferroni' and args.multipleTesting != 'sidak' and args.multipleTesting != 'holm-sidak' and args.multipleTesting != 'holm' and args.multipleTesting != 'simes-hochbergh' and args.multipleTesting != 'hommel' and args.multipleTesting != 'fdr_bh' and args.multipleTesting != 'fdr_by' and args.multipleTesting != 'fdr_tsbh' and args.multipleTesting != 'fdr_tsbky':
            sys.exit('ERROR! input for multipleTesting cannot be recognized, please check input again')
        if float(args.criticalValue) <= 0.0 or float(args.criticalValue) > 1.0:
            sys.exit('ERROR! input for criticalValue cannot be recognized, please check input again')
        test_names = args.condition1
        test_names = test_names.split(',')
        control_names = args.condition2
        control_names = control_names.split(',')
        if len(test_names) < 2 or len(control_names) < 2:
            sys.exit("ERROR! Please include at least two samples for each conditions")

        # extract names and contruct a list of list with samples names and conditions or types
        for i in test_names:
            sample_types += [[i, 'test']]

        for i in control_names:
            sample_types += [[i, 'control']]

        stat_analysis_main_function(sample_types, args.output, args.biosignificant, args.minimalGeneCount, args.minimalExonCount, args.criticalValue, args.outlierDetection, args.outlierMethod, args.multipleTesting)
