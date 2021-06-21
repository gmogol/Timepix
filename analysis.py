import numpy as np
import matplotlib.pyplot as plt

def to_plot_2dhist(signal_df,var1,var2,*,full_output=False, **kwargs):
    '''
    Plots a 2d histogram of timepix data
    Args:
        var1 : string
            first variable to plot. allowed values are x,y and delta_t
        var2 : String
            second variable to plot. allowed values are x,y and delta_t
        full_output : bool, optional
            if True returns bin edges as well as binned counts
    Kwargs:
        bin1 : int or array_like, optional
            bins for the first variable if nothing is given
            default bins are chosen for the variable

        bin2 : int or array_like , optional
            bins for the second variable if nothing is given
            default bins are chosen for the variable
        signal_factor : float, optional
            factor to multiply the signal counts
        background_factor : float, optional
            factor to multiply the background counts

            
    Returns:
        bin_counts : array_like
            directly plottable bin_counts using plt.imshow
        bin_edges : array_like, optional
            bin edges, only returned if full_output is True
    '''

    background = kwargs.get('background')
    bin1 = kwargs.get('bin1')
    bin2 = kwargs.get('bin2')
    bint = kwargs.get('bint')
    signal_factor = kwargs.get('signal_factor',1)
    background_factor = kwargs.get('background_factor',1)
    

    


        
    if var1 not in ['x','y','delta_t'] or var2 not in ['x','y','delta_t']:
        raise Exception('Only allowed variables are x,y and delta_t')
        
    def choose_bins(var):
        if var in ['x','y']:
            bins = np.linspace(0,255,100)
        if var == 'delta_t':
            if bint is not None:
                bins = bint
            else:
                bins = np.linspace(4000,5000,100,endpoint=False)
        return bins
    
    if bin1 is None:
        bin1 = choose_bins(var1)
    if bin2 is None:
        bin2 = choose_bins(var2)
        
    signal_hist = np.histogram2d(signal_df[var1],
                                 signal_df[var2],
                                 bins=(bin1,bin2))
    if background is not None:
        background_hist = np.histogram2d(background[var1],
                                     background[var2],
                                     bins=(bin1,bin2))
    if full_output:
        print(full_output)
        return (signal_factor*signal_hist[0]-background_factor*background_hist[0]), signal_hist[1:]
    if background is None:
        return signal_hist[0].T
    return signal_factor*signal_hist[0].T-background_factor*background_hist[0].T


def plot_integrated(signal,*,cmap='viridis',**kwargs):

    plot_file = kwargs.get('plot_file')
    title = kwargs.get('title')
    
    fig = plt.figure(figsize=(15,10))
    if title is not None:
        fig.suptitle(title)
    fig.add_subplot(231)
    plt.imshow(to_plot_2dhist(signal,'x','delta_t',**kwargs)
               ,aspect='auto',cmap=cmap)
    plt.xlabel('x')
    plt.ylabel('time')



    fig.add_subplot(233)
    plt.imshow(to_plot_2dhist(signal,'y','delta_t',**kwargs)
               ,aspect='auto',cmap=cmap)
    plt.xlabel('y')
    plt.ylabel('time')

    fig.add_subplot(234)
    plt.plot(to_plot_2dhist(signal,'x','delta_t',**kwargs).sum(axis=0))
    plt.xlabel('x')
    plt.ylabel('counts')

    fig.add_subplot(235)
    plt.xlabel('time')
    plt.plot(to_plot_2dhist(signal,'y','delta_t',**kwargs).sum(axis=1))

    fig.add_subplot(236)
    plt.xlabel('y')
    plt.plot(to_plot_2dhist(signal,'y','delta_t',**kwargs).sum(axis=0))

    if plot_file is not None:
        plt.savefig(plot_file)

    plt.show()


def plot_projections(signal,*,cmap='viridis',**kwargs):

    title = kwargs.get('title')
    plot_file = kwargs.get('plot_file')



    fig = plt.figure(figsize=(10,10))
    if title is not None:
        fig.suptitle(title)
    
    fig.add_subplot(221)




    plt.imshow(to_plot_2dhist(signal,'x','y',**kwargs),
               aspect='auto',
               cmap=cmap)
    plt.xlabel('x')
    plt.ylabel('y')


    fig.add_subplot(223)



    plt.imshow(to_plot_2dhist(signal,'x','delta_t',**kwargs),
               aspect='auto',
               cmap = cmap)
    plt.xlabel('x')
    plt.ylabel('time')


    fig.add_subplot(222)


    plt.imshow(to_plot_2dhist(signal,'delta_t','y',**kwargs),
               aspect='auto',
               cmap = cmap)
    plt.xlabel('time')
    plt.ylabel('y')
    if plot_file is not None:
        plt.savefig(plot_file)
    plt.show()
