import numpy as np
import random

def make_isodomain_map(nx, ny, gx=771, gy=771):
    return np.ones((nx-1, ny)) * gx, np.ones((nx, ny-1)) * gy

def make_isodomain_map_heter(nx, ny, randcnx=0.1, randcny=0.1, gx=700, gy=700, gxcn=7000, gycn=7000):
    weightsx, weightsy = np.ones((nx, ny-1)) * gx, np.ones((nx-1, ny)) * gy
    for i in range(len(weightsx)):
        for j in range(len(weightsx[0])):
            x_rand = np.random.random(size=())
            if x_rand <= 0.1:
                weightsx[i][j] = gxcn
    for i in range(len(weightsy)):
        for j in range(len(weightsy[0])):
            y_rand = np.random.random(size=())
            if y_rand <= 0.1:
                weightsy[i][j] = gycn
    return weightsx, weightsy


def make_isodomain_map_heter_cut(nx, ny, randcnx=0.1, randcny=0.1, gx=700, gy=700, gxcn=7000, gycn=7000, g_cut=0):
    weightsx, weightsy = np.ones((nx, ny-1)) * gx, np.ones((nx-1, ny)) * gy
    for i in range(len(weightsx)):
        for j in range(len(weightsx[0])):
            x_rand = np.random.random(size=())
            if x_rand <= 0.1:
                weightsx[i][j] = gxcn
            if i <= len(weightsx)//2 and len(weightsx[0])//2 <= j <= len(weightsx[0])//2 + 2:
                weightsx[i][j] = g_cut
    for i in range(len(weightsy)):
        for j in range(len(weightsy[0])):
            y_rand = np.random.random(size=())
            if y_rand <= 0.1:
                weightsy[i][j] = g_cut
            if i <= len(weightsx)//2 and len(weightsx[0])//2 <= j <= len(weightsx[0])//2 + 2:
                weightsy[i][j] = g_cut

    # for i in range(ny//2, ny-2):
    #     for j in range(nx//2, nx//2+10):
    #         print(i,j)
    #         weightsy[i][j] = 0
    #         weightsx[i][j] = 0

    return weightsx, weightsy


def crop_center(img, startx, starty, sizex, sizey):
    # y,x = img.shape
    # startx = x//2 - cropx//2
    # starty = y//2 - cropy//2    
    return img[int(starty):int(starty+sizey), int(startx):int(startx+sizex)]


def find_bounds(ctags):
    x_bounds = (np.diff(ctags, axis=1) == 0).astype('uint8')
    y_bounds = (np.diff(ctags, axis=0) == 0).astype('uint8')
    
    return x_bounds, y_bounds


def find_contacts(contacts):
    x_contacts = (np.diff(contacts, axis=1) != 0).astype('uint8')
    y_contacts = (np.diff(contacts, axis=0) != 0).astype('uint8')
    
    return x_contacts, y_contacts

def find_contacts_adaptive(contacts, ctags, ratio=0.5):
    x_bounds = (np.diff(ctags, axis=1) == 0).astype('uint8')
    y_bounds = (np.diff(ctags, axis=0) == 0).astype('uint8')

    x_contacts = (np.diff(contacts, axis=1) != 0).astype('uint8')
    y_contacts = (np.diff(contacts, axis=0) != 0).astype('uint8')
    
    return x_contacts, y_contacts


def find_cms(ctags, types):
    x_cms = np.zeros((ctags.shape[0], ctags.shape[1] - 1))
    y_cms = np.zeros((ctags.shape[0] - 1, ctags.shape[1]))
    
    for i in range(ctags.shape[0]):
        for j in range(ctags.shape[1] - 1):
            if types[ctags[i, j+1]] - types[ctags[i, j]] == 0 and types[ctags[i, j]] == 1:
                x_cms[i, j] = 1
    
    for i in range(ctags.shape[0] - 1):
        for j in range(ctags.shape[1]):
            if types[ctags[i+1, j]] - types[ctags[i, j]] == 0 and types[ctags[i, j]] == 1:
                y_cms[i, j] = 1
                
    return x_cms, y_cms


def make_potts_map(ctags_path='VCT/output/ctags3000.out', conts_path='VCT/output/contactM3000.out', types_path='VCT/output/types.out',
                  g_in=771, g_gap=500, g_bound=10, startx=1, starty=1, sizex=10, sizey=10):
    ctags = crop_center(np.loadtxt(ctags_path), startx, starty, sizex, sizey).astype('uint32')
    conts = crop_center(np.loadtxt(conts_path), startx, starty, sizex, sizey).astype('uint32')
    types = np.loadtxt(types_path)
    print(types)
                   
    bound_mask_x = (np.diff(ctags, axis=1) != 0).astype('uint8')
    bound_mask_y = (np.diff(ctags, axis=0) != 0).astype('uint8')
    
    
    x_cm, y_cm = find_cms(ctags, types)
    x_b, y_b = find_bounds(ctags)
    x_cont, y_cont = find_contacts(conts)
    
    gx = x_b * g_in * x_cm + g_gap * x_cont * (x_b == 0).astype('uint8') + (x_b == 0).astype('uint8') * g_bound
    gy = y_b * g_in * y_cm + g_gap * y_cont * (y_b == 0).astype('uint8') + (y_b == 0).astype('uint8') * g_bound

    return gx, gy   
