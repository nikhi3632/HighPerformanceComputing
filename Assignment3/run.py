#!/usr/bin/env python3

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import PIL.Image as pimg

def plot_dat_file_and_save(dat_file):
    cwd = os.getcwd()
    file_path = cwd + '/' + dat_file
    x, h_a, h_exa = np.loadtxt(file_path, unpack=True)
    fig, ax = plt.subplots(figsize=(10,5))
    ax.plot(x, h_a, 'ro', label= 'Numerical')
    ax.plot(x, h_exa, 'g*', label= 'Exact')
    ax.legend()
    ax.legend(loc=2)
    plt.xlim([-8,8])
    plt.ylim([-0.4,1.2])
    image_name = dat_file[:-4] + '.png'
    fig.savefig(image_name)
    image_path = cwd + '/' + image_name
    image = pimg.open(image_path)
    im = image.convert('RGB')
    im.save(dat_file[:-4] + '.pdf')

def main():
    try:
        file_name = 'shallow.dat'
        plot_dat_file_and_save(file_name)
    except Exception as e:
        print('Exception Occured: ', e)
        sys.exit(1)

if __name__ == '__main__':
    main()