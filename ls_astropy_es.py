# -*- coding: utf-8 -*-
"""
Plot graph according to the DAT file 
@author: Danil Borchevkin
"""

import csv
import glob
import os
import matplotlib.pyplot as plt
import scipy.signal as signal
import numpy as np
import math
from astropy.timeseries import LombScargle

def read_raw_file_data(filepath):
    '''
    Read data in list
    '''

    raw_data = list()

    # Get raw data
    with open(filepath, 'r') as dest_f:
        raw_data = dest_f.readlines()

    return raw_data

def process_file(data, out_filepath, window, step):
    line_cursor = 0

    while (line_cursor < (len(data) - window)):
        with open(out_filepath + '_c' + "{:08d}".format(line_cursor) + '_w' + str(window) + '_s' + str(step) + ".dat", 'w') as dest_f:
            for i in range(window):
                dest_f.write(data[line_cursor + i])
        line_cursor += step

def read_file_data(filepath):
    '''
    Read data in [[val,time],[val, time]] format
    '''

    raw_data = None
    data = list()

    # Get raw data
    with open(filepath, 'r') as dest_f:
        data_iter = csv.reader(dest_f,delimiter="\t")
        raw_data = [raw_data for raw_data in data_iter]

    # Convert data to list. If data is absent set it to None
    for raw_val in raw_data:
        amp = 0
        time = 0
        try:
            amp = float(raw_val[0])
        except:
            amp = None
        finally:
            time = float(raw_val[1])
            data.append([amp, time])

    return data

def save_to_ascii_file(data_list, out_filepath, header=[]):
    '''
    Save data in format [[],[]] into DAT file 
    - CSV 
    - with \t delimeter 
    - \n line endings
    '''
    write_list = []

    for data in data_list:
        output_str = ""
        for val in data:
            output_str += str(val) + "\t"
        output_str = output_str[:-1]
        output_str += "\n"
        write_list.append(output_str)

    with open(out_filepath,"w") as f:
        f.writelines(write_list)

def plot_graph(data, out_filepath, to_display=False, save_to_disk=True):
    '''
    Plot grapth and return its data

    Params
    data - input data in list of lists with pair value and time
    out_filepath - out file name path for create
    to_display - if set to true then graph will be shown on the display
    save_to_disk - if set to true then graph will be saved on the disk

    Return
    List of lists of graph values in form [freq, period, pgram_value, time_value]
    '''

    output_data = list()

    x = list()
    y = list()

    # Get first time value as constant time value for all window
    time_value = data[0][1]

    for val_pair in data:
        if val_pair[0] != None:
            x.append(val_pair[1])
            y.append(val_pair[0])

    # Calculate Lomb-Scargle periodogram Astropy
    astropy_pgram = LombScargle(x, y, normalization='psd')
    astropy_freq, astropy_power = astropy_pgram.autopower()
    astropy_false_alarm_probability = astropy_pgram.false_alarm_probability(astropy_power.max(), method='baluev')

    # Create figure with 2 subplots
    fig = plt.figure()
    source_ax = fig.add_subplot(211)
    astropy_pgram_ax = fig.add_subplot(212)

    #Now make a plot of the input data:
    source_ax.plot(x, y, 'b+')

    # astropy periodogram
    astropy_pgram_ax.plot(astropy_freq, astropy_power,'g')
    astropy_pgram_ax.text(0.95, 0.95, "FAP(first_peak) = {:.4f}%".format(astropy_false_alarm_probability),
        verticalalignment='top', horizontalalignment='right',
        transform=astropy_pgram_ax.transAxes,
        color='green', fontsize=15)

    if to_display:
        plt.show()

    if save_to_disk:
        plt.savefig(out_filepath)

    # Generate output
    for idx, freq in enumerate(astropy_freq):
        period = 1 / freq
        output_data.append([freq, period, astropy_power[idx], time_value])

    plt.cla()
    plt.clf()
    plt.close(fig)

    return output_data

def process_windowed_files(path, output_file_path):
    files = glob.glob(path + "*.dat")  

    for filepath in files:
        # Reject old merged files
        if "!" in filepath:
            continue

        # Reject old windowed files
        if "windowed" in filepath:
            continue

        print("Process >> " + filepath)


        read_data = read_file_data(filepath)
        out_dat_filepath = path + os.path.basename(filepath) + "_windowed" + ".dat"
        out_png_filepath = path + os.path.basename(filepath) + "_windowed" + ".png"

        output_data = plot_graph(read_data, 
                                out_png_filepath)

        print("Saved PNG to >> " + out_png_filepath)
        save_to_ascii_file(output_data, out_dat_filepath)
        print("Saved DAT to >> " + out_dat_filepath)

        try:
            os.remove(output_file_path)
        except Exception as e:
            pass
        finally:
            pass

        windowed_files = glob.glob(path + "*_windowed.dat")
        for windowed_file in windowed_files:           
            with open(windowed_file, 'r') as windowed_f:
                data = windowed_f.read()
                with open(output_file_path, 'a') as merged_file:
                    merged_file.write(data)

def main():
    print("Script is started")

    files = glob.glob("./input/*.dat")              # Change path here or write filepath
    OUTPUT_PATH = "./output/"                       # Change output here
    WINDOW = 648                                    # Change window value here
    STEP = 24                                      # Change step value here
    FREQ_START = 0.08                               # Change freq start here
    FREQ_END = 1.0                                  # Change freq end here
    FREQ_NUM = 500                               # Change freq num here

    for filepath in files:
        print("Process >> " + filepath)

        read_lines = read_raw_file_data(filepath)
        out_dat_filepath = OUTPUT_PATH + os.path.basename(filepath)
        process_file(read_lines, out_dat_filepath, WINDOW, STEP)
        process_windowed_files(OUTPUT_PATH, f'{OUTPUT_PATH}!{os.path.basename(filepath)}_merged_file.dat')

        print(f"<{filepath}> succesful processed by the script")

    print("Script is finished")

if __name__ == "__main__":
    main()