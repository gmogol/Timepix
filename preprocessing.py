import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO

def tdc_stripper(read_filename, write_filename,total_lines = -1):
    """
    Args:
    read_filename : string
            the file to strip off the TDC data
    write_filename : string
	    the file to write the TDC data
	total_lines : int (optional)
	    total number of lines to read, default = -1 = inf
    Returns:
	None
    """

    line_counter = 0
    with open(read_filename,'r') as rfile:
        with open(write_filename,'w') as wfile:
            for line in rfile:
                line_counter+=1
                if line_counter == total_lines:
                    break
                if line[:3] == 'TDC':
                    wfile.write(line)


def preproces_tdcs(filename=None,debug=True,tdc_df=None):
    ''' 
    Args:
        filename : string
            name of the file
        debug : bool
              if True returns all dataframes
        tdc_df : pd.DataFrame Object
            instead of reading from file one can give the df directly
    Returns:
        time_df : pd.DataFrame Object
            
    It takes the data file for TDC's ONLY, which is assumed to be in the form
    TDCn: t0.
    We first seperate TDC's apart and then merge them, which creates nulls.
    The nulls for TDC2 are back filled i.e. if a cell N is null then it's filled 
    with the last cell <N which isn't null. TDC1 nulls aren't filled but rather dropped
    because we are really interested in TDC2. Finally we calculate the time difference
    between TD1 and TDC2. Finally we gate the difference by 5*10^4 ns = 50 us.
    '''
    if filename is not None:
        tpx_df =pd.read_csv(filename,
                    sep=':',
                    error_bad_lines=False,
                    warn_bad_lines=False,
                    header=None)
    elif tdc_df is not None:
        tpx_df = tdc_df
    else:
        raise(Exception('Either a filename or a tdc data frame needs to be given!'))
    tpx_df.columns = ['tdc','time']

    print('tpx_df nulls\n',tpx_df.isna().sum())
    tdc_one = tpx_df[tpx_df.tdc == 'TDC1']
    tdc_two = tpx_df[tpx_df.tdc == 'TDC2']
    time_df = tdc_one[['time']].join(tdc_two.time,how='outer',lsuffix='1',rsuffix='2')
    time_df['time2'] = time_df['time2'].fillna(method='bfill')
    time_df['time1'] = time_df['time1'].fillna(method='bfill')

    #time_df.dropna(inplace=True)
    time_df = time_df.astype(np.float64)
    time_df['delta_t'] = (time_df.time2-time_df.time1)*1e9
    time_df.reset_index(inplace=True)
    time_df.drop('index',axis=1,inplace=True)
    #time_df = time_df[np.abs(time_df.delta_t)< 5e3]
    if debug:
        return tpx_df, tdc_one, tdc_two, time_df
    else:
        print(time_df.head(10))

        return time_df


def final_processing(df_str,debug=False):
    if debug:
        tpx_data = df_str
    else:
        tpx_data = pd.read_csv(StringIO(df_str),index_col=False)
    tpx_data['totx'] = tpx_data['tot']*tpx_data['x']
    tpx_data['toty'] = tpx_data['tot']*tpx_data['y']
    tpx_data['delta_t'] = (tpx_data['tdc2']-tpx_data['tdc1'])
    tpx_data['ones'] =1
    batches = tpx_data.groupby('batch').sum()
    batches['delta_t'] = batches['delta_t']/batches['ones']
    batches['x'] = batches['totx']/batches['tot']
    batches['y'] = batches['toty']/batches['tot']
    batches.drop(['toa','tot','tdc1','tdc2','totx','toty','ones'],axis=1,inplace=True)
    if debug:
        return batches
    return batches.to_csv(header=False)




    

def preprocess(data_file,tdc_file,write_file, **kwargs):
    file_check = input('You will write on the file ' + write_file + ' continue? y/n')
    if file_check == 'n':
        return
    
    _,tdc_one, tdc_two,_= preproces_tdcs(tdc_file,debug=True)
    tdc_one = np.array(tdc_one['time'])
    tdc_two = np.array(tdc_two['time'])


    NTOTBATCHES = kwargs.get('ntot_batches',-1)    # Max number of batches to process -1=inf
    NBATCH = kwargs.get('nbatch',2_000)    # Number of batches to collect before writing it to the file
    BEGINLINE= kwargs.get('beginline',0)   # starting line, used mainly for debugging
    line_counter = 0
    NLINES = kwargs.get('nlines',-1)  # Max number of lines to be read -1=inf
    batch = 0

    tdc_updated = False

    t1_eps = kwargs.get('t1_eps',9e-4) #tdc1 tolerance [s]
    t2_eps = kwargs.get('t2_eps',3e-6) #tdc2 tolerance [s]

    x_threshold = kwargs.get('x_threshold',15)  # max batch size in x direction for rejecting many electrons
    y_threshold = kwargs.get('y_threshold',15) # max batch size in y direction for rejecting many electrons
    
    xmax = 0
    xmin = 1000 # detector has only 256x256 pixels
    ymax = 0
    ymin = 1000 # detector has only 256x256 pixels

    batch_lines = ''
    df_str = ''

    many_electron_reject = 0
    tdc1_reject = 0
    tdc2_reject= 0

    # We gate tdc2-tdc1 differnce (which corresponds to cable length)
    # and reject if this difference is too much or too little
    # physical_max/min_t's were determined by an preliminary analysis of the data
    # and they are given in ns.

    unphysical_reject = 0

    gate_time = kwargs.get('gate_time',False)
    physical_max_t = kwargs.get('physical_max_t',0)   # max allowed tdc2-tdc1 difference [ns]
    physical_min_t = kwargs.get('physical_min_t',0)   # min allowed tdc2-tdc1 difference [ns]

    toa_overflow_counter = 0
    tdc1_overflow_counter = 0
    tdc2_overflow_counter = 0

    toa_overflow_threshold = (25*2**30/1e9)   # when pixel's toa overflows 
    tdc_overflow_threshold = (25*2**32/1e9)  # when tdcs overflow

    tdc1_overflow = False
    tdc2_overflow = False
    toa_overflow = False

    overflow_threshold = 10     # sets after which we start checking overflow
                                # i.e. difference bw subsequent values need to
                                # be at least this threshold

    toa = 0



    tdc1_iter = tdc_one.flat   # flat produces an iterator from np.array
    tdc1 = next(tdc1_iter)
    tdc1_next = next(tdc1_iter)
    tdc1_batch = 0 # chosen tdc1 for the batch


    tdc2_iter = tdc_two.flat   # flat produces an iterator from np.array
    tdc2 = next(tdc2_iter)
    tdc2_next = next(tdc2_iter)
    tdc2_batch = 0 # chosen tdc2 for the batch


    header = 'batch,x,y,delta_t\n'
    counter_list = ['toa_overflow_counter' ,
                    'tdc1_overflow_counter' ,
                    'tdc2_overflow_counter' ,
                    'unphysical_reject'  ,
                    'many_electron_reject',
                    'tdc1_reject'  ,
                    'tdc2_reject',  
                    'batch'  ]



    with open(write_file,'w') as wfile:
        wfile.write('')

    import re
    re_search = re.compile('\w+:')
    with open(data_file,'r') as rfile:
        with open(write_file,'a') as afile:
            afile.write(header)

            def overflow_corr(num,num_counter,num_threshold): #corrects for overflow
                return num + num_counter*num_threshold

            for line in rfile:
                line_counter+=1
                if line_counter < BEGINLINE:
        #            toa_overflow_counter =3 #for 12_000_000
                    continue

                if line_counter %1_000_000 == 0:
                    print(f'{line_counter:_} lines read')


                if batch % NBATCH == 0 and batch>0 and df_str != '':
                    df_str='toa,tot,x,y,tdc1,tdc2,batch\n'+df_str
                    to_write = final_processing(df_str)
                    afile.write(to_write)
                    notes_file_name = write_file.split('.')[0] + '_notes'
                    with open(notes_file_name,'w') as notes:
                        for counter in counter_list:
                            notes.write(counter + ' : ' + str(eval(counter)) + '\n')
                    to_write = ''
                    #break
                    df_str = ''
                    print(f'\nBatch numer {batch:_} is written')
                    print(f'{line_counter:_} lines read')
                    print('-----')



                if batch == NTOTBATCHES:
                    print('NBATCHES reached')
                    break
                if line_counter == NLINES:
                    print('NLINES reached')
                    break


    ############################################ TOA #######################################


                str_list = line.split('TOA:')
                try:
                    toa_prev = toa
                    toa = np.float64(str_list[1].split(',')[0])
                    toa = overflow_corr(toa,toa_overflow_counter,toa_overflow_threshold)

                    if toa_prev > toa + overflow_threshold:

                        toa_overflow_counter += 1
                        #print('toa_prev %s toa %s toa_prev-toa %s' % (toa_prev,toa,toa_prev-toa))
                        toa += toa_overflow_threshold
                        print('### TOA Overflow NO %s' %toa_overflow_counter )
                        print('toa_prev %s toa %s tdc1 %s tdc2 %s'  % (toa_prev,toa,tdc1,tdc2))



                except IndexError: # This only happens for tdc lines
                    continue


    ############################################ TDC2  #######################################


                # tdc2 iterator is only updated when the difference of toa with the current tdc2
                # is greater than the difference of toa with the next tdc2
                while abs(toa-tdc2) >= abs(toa-tdc2_next):
                    tdc2 = tdc2_next
                    try:
                        tdc2_next = next(tdc2_iter)
                    except StopIteration:
                        print('\n #########______######### \n')
                        for counter in counter_list:
                            print(f'{counter}: {eval(counter)}')
                        return None

                    


                    # Handle overflow
                    tdc2_next = overflow_corr(tdc2_next,tdc2_overflow_counter,tdc_overflow_threshold)

                    if tdc2 > tdc2_next+overflow_threshold:
                       # print('tdc2 %s tdc2_next %s tdc1 %s toa %s' % (tdc2,tdc2_next ,tdc1,toa))

                        tdc2_overflow_counter += 1 
                        tdc2_next += tdc_overflow_threshold
                        print('### TDC2 Overflow NO %s' %tdc2_overflow_counter )
                        print('tdc2 %s tdc2_next %s tdc1 %s toa %s' % (tdc2,tdc2_next ,tdc1,toa))





                    tdc_updated = True


    ############################################ MANY ELECTRON  #######################################

                # We want to gate >1 electron events. This is achieved by putting a max-min threshold
                # on x and y. 

                # before we deal with many electrons we want to reject unphysical events. these are
                # due to e.g. dark counts and faults of the detector

                if tdc_updated:
                    # *1e9 is to convert from s to ns
                    if ((tdc2_batch-tdc1_batch)*1e9 < physical_max_t and \
                        (tdc2_batch-tdc1_batch)*1e9 > physical_min_t) or (not gate_time):
                        if abs(xmax-xmin) < x_threshold and abs(ymax-ymin) < y_threshold:
                            df_str = df_str + batch_lines
                            batch += 1
                        else:
                            if np.random.rand()<1/1250:   # arbitrary number that determines when to plot pics
                                many_e = pd.read_csv(StringIO('toa,tot,x,y,tdc1,tdc2,batch\n'+batch_lines))
                                try:
                                    many_e.plot(x='x',y='y',linestyle='None',marker='o')
                                    plt.show()
                                except:
                                    None
                            many_electron_reject +=1
                            if many_electron_reject % 5_000==0:
                                print(f'{many_electron_reject:_} physical batches were \
                                      rejected due to many electrons')
                    else:
                        unphysical_reject += 1
                        if unphysical_reject % 5_000 == 0:
                            print(f'{unphysical_reject:_} batches were rejected due time gating')



                    tdc_updated = False
                    batch_lines =''
                    xmax = 0
                    xmin = 1000
                    ymax = 0
                    ymin = 1000

    ############################################ TDC1  #######################################


                # when the difference bw tdc1 and toa exceeds the threshold we need to update tdc1
                while toa-tdc1>=t1_eps:
                    tdc1 = tdc1_next
                    try:
                        tdc1_next = next(tdc1_iter)
                    except StopIteration:
                        print('\n #########______######### \n')
                        for counter in counter_list:
                            print(f'{counter}: {eval(counter)}')
                        return None
                    # Handle overflow
                    tdc1_next = overflow_corr(tdc1_next,tdc1_overflow_counter,tdc_overflow_threshold)

                    if tdc1 > tdc1_next + overflow_threshold:
                        tdc1_overflow_counter += 1
                       # print('tdc1 %s tdc1_next %s tdc2 %s toa %s' % (tdc1,tdc1_next ,tdc2,toa))

                        tdc1_next += tdc_overflow_threshold
                        print('### TDC1 Overflow NO %s' %tdc1_overflow_counter )
                        print('tdc1 %s tdc1_next %s tdc2 %s toa %s' % (tdc1,tdc1_next ,tdc2,toa))




                if toa-tdc1 >=0 and toa-tdc1<t1_eps:  
                    # We assume that tdc1 signal must be the first signal hence >0
                    # if the difference between time of arrival (toa) and tdc1 signal
                    # is below threshold and toa-tdc2 is below threshold we write into the file
                    if abs(toa-tdc2) < t2_eps:
                        line = ''.join(re_search.split(line))[:-1] + (', %s , %s , %s\n' % (tdc1,tdc2,batch))
                        tdc1_batch = tdc1
                        tdc2_batch = tdc2
                        batch_lines = batch_lines + line
                        x,y = line.split(',')[2:4]
                        x,y = int(x),int(y)
                        xmax = max(x,xmax)
                        xmin = min(x,xmin)
                        ymax = max(y,ymax)
                        ymin = min(y,ymin)
                    else:
                        tdc2_reject +=1
                        if tdc2_reject %100_000 == 0:
                            print(f'{tdc2_reject:_} lines were rejected due to tdc2 mismatch')
                else:
                    tdc1_reject += 1
                    if tdc1_reject %100_000 == 0:
                            print(f'{tdc1_reject:_} lines were rejected due to tdc1 mismatch')

    print('\n #########______######### \n')
    for counter in counter_list:
        print(f'{counter}: {eval(counter)}')
    return None
