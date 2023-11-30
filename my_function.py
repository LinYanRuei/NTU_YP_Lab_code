import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import os
from PIL import Image
import pwlf

def fit_func(x,a,b):
    return a*x+b


def calculate_Bandgap(number:int,list,list_sample_bias):

    """
    This is a function to get CBM and VBM
    number: order of title
    list: list for di/dv value 
    list_sample_bias: list for bias voltage
    return is CBM & VBM

    """

    offset_didv = list
    self_def = 40

    slope_val = np.inf
    k = 0
    didv_tran = []
    for i in range (len(offset_didv)):
        if i < self_def: continue
        slope = offset_didv[i-self_def]-offset_didv[i]
        if np.abs(slope) < slope_val: 
            slope_val = np.abs(slope)
            k = i-self_def
            didv_tran = offset_didv[(k-self_def):k]
    std = np.std(didv_tran)
    mean = np.mean(didv_tran)


    didv_left = []
    didv_right = []
    bias_left = []
    bias_right = []
    for i in range(len(offset_didv)):
        if i <= k and offset_didv[i] >= mean+2*std:
            didv_left.append(offset_didv[i])
            bias_left.append(list_sample_bias[0][i])
        if i >= k + self_def and offset_didv[i] >= mean+2*std:
            didv_right.append(offset_didv[i])
            bias_right.append(list_sample_bias[0][i])

    popt_l, pcov_l = fit(fit_func, bias_left, didv_left,maxfev = 1000000)
    popt_r, pcov_r = fit(fit_func, bias_right, didv_right,maxfev = 1000000)
    bias_left = np.array(bias_left)
    bias_right = np.array(bias_right)

    residuals = didv_left - fit_func(bias_left, *popt_l)
    ss_residuals = np.sum(residuals**2)
    ss_total = np.sum((didv_left - np.mean(didv_left))**2)
    r_squared_l = 1 - (ss_residuals / ss_total)
    

    residuals = didv_right - fit_func(bias_right, *popt_r)
    ss_residuals = np.sum(residuals**2)
    ss_total = np.sum((didv_right - np.mean(didv_right))**2)
    r_squared_r = 1 - (ss_residuals / ss_total)
    
    for i in range(len(didv_left)):
        if r_squared_l >= 0.93:
            continue
        else:
            if i%7 == 0:
                bias_left = bias_left.tolist()
                bias_left.pop()
                didv_left.pop()
            else:
                bias_left = bias_left.tolist()
                didv_left.pop(0)
                bias_left.pop(0)
            popt_l, pcov_l = fit(fit_func, bias_left, didv_left,maxfev = 1000000)
            bias_left = np.array(bias_left)

            residuals = didv_left - fit_func(bias_left, *popt_l)
            ss_residuals = np.sum(residuals**2)
            ss_total = np.sum((didv_left - np.mean(didv_left))**2)
            r_squared_l = 1 - (ss_residuals / ss_total)

    for i in range(len(didv_right)):
        if r_squared_r >= 0.97:
            continue
        else:
            if i%3 == 0:
                bias_right = bias_right.tolist()
                didv_right.pop()
                bias_right.pop()
            else:
                bias_right = bias_right.tolist()
                bias_right.pop(0)
                didv_right.pop(0)
            popt_r, pcov_r = fit(fit_func, bias_right, didv_right,maxfev = 1000000)
            bias_right = np.array(bias_right)

            residuals = didv_right - fit_func(bias_right, *popt_r)
            ss_residuals = np.sum(residuals**2)
            ss_total = np.sum((didv_right - np.mean(didv_right))**2)
            r_squared_r = 1 - (ss_residuals / ss_total)


    print('VBM = ',(offset_didv[k]-popt_l[1])/popt_l[0])
    print('CBM = ',(offset_didv[k]-popt_r[1])/popt_r[0])
    print('Bandgap = ',((offset_didv[k]-popt_r[1])/popt_r[0]-(offset_didv[k]-popt_l[1])/popt_l[0]))
    print(r_squared_l,r_squared_r)

    plt.plot(list_sample_bias[0][k+self_def:],fit_func(list_sample_bias[0][k+self_def:],*popt_r),'r')
    plt.plot(list_sample_bias[0][:k],fit_func(list_sample_bias[0][:k],*popt_l),'r')
    plt.axhline(y=offset_didv[k],c='black',linestyle='dashed')
    plt.scatter(list_sample_bias[0],offset_didv,s=8)
    plt.tick_params(direction='in')
    plt.ylabel('di/dv')
    plt.xlabel('Bias(V)')
    plt.title(f'CBMVBM_{number}')
    
    # plt.savefig(f'CBMVBM_{number}.png')
    plt.show()

    return (offset_didv[k]-popt_r[1])/popt_r[0],(offset_didv[k]-popt_l[1])/popt_l[0]


def open_csv_file(path: str = '', name: str = '', blank_line: int = 13):
    """
    Open CSV file

    path: data file path
    name: data name(e.g. XXX.csv), must include .csv
    blank_line: top line that should be removed, default is 13
    return: 2D List
    """


    if path != '': 
        full_name = os.path.join(path, name)
    else: 
        full_name = name
    with open(full_name, newline='',encoding="utf-8") as csvfile:  ## di/dv data
        csv_reader = csv.reader(csvfile)
        list_reporter = list(csv_reader)
        list_reporter = np.array(list_reporter).T 
        list_reporter = list_reporter[:,blank_line:]
        list_result = np.rot90(list_reporter,1).astype(float)
    return list_result


def height_linecut(y_value: int, x_start: int = 5, x_end: int = 495, x_length= int, measure_array= None):
    avg_list = []
    y_start = y_value   
    x_range = x_end-x_start
    label = x_end
    result = np.array(measure_array)
    i = np.linspace(x_start,x_end,x_range)
    plt.plot(i,result[y_start][x_start:x_end])
    plt.show()
    
    
    fig, ax = plt.subplots()
    plt.plot(i,result[y_start][x_start:x_end])
    x_end = 0
    print('close at input = ', label)
    while x_end < label :
        x_start = int(input('left of truncate = '))
        x_end = int(input('right of truncate = '))
        if x_start == 0:
            continue
        x_range = x_end- x_start
        g, m = 0,0
        for n in range(x_range):
            g = result[y_start,n+x_start]
            m = m+float(g)
        avg = m/x_range
        avg_list.append(avg)
        i = np.linspace(x_start,x_end,x_range)
        result = result.astype(float)
        j = avg*np.ones(len(i))
        
        plt.plot(i,j,'r')
    
    xtick_positions = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 10)
    ytick_positions = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 10)
    ax.set_xticks(xtick_positions)
    ax.set_yticks(ytick_positions)
    # Set the tick labels for the x-axis and y-axis
    ax.set_xticklabels(np.linspace(0, x_length, 10).astype(int))

    ax.tick_params(axis='both', direction='in', length=6, width=1, colors='k')
    plt.ylabel('y (nm)')
    plt.xlabel('x (nm)')
    plt.show()
    return(avg_list)



def RGB2Gray(array):
    image = Image.fromarray(array, 'RGB')

    # 轉換為灰階
    gray_image = image.convert('L')

    # 如果需要，將灰階圖像轉換回 numpy array
    gray_array = np.array(gray_image)
    return gray_array

def STS_Fitting(x,y):
    """
    This is a function to calculate STS bandgap
    X: STS sample bias
    Y: STS data
    return: Bandgap (float)
    """
    my_pwlf = pwlf.PiecewiseLinFit(x, y)
    breaks = my_pwlf.fit(3)  # 3表示三段線性
    # x0 = np.array([min(x), -1.45, 0.73, max(x)])
    # breaks = my_pwlf.fit_with_breaks(x0)
    x_hat = np.linspace(min(x), max(x), 100)
    y_hat = my_pwlf.predict(x_hat)

    # 繪圖
    plt.scatter(x, y, label='Data')
    plt.plot(x_hat, y_hat, color='black', label='Piecewise Linear Fit')


    for breakpoint in breaks[1:-1]:  # 排除起始和結束點
        plt.axvline(x=breakpoint, color='red', linestyle='--')
        plt.text(breakpoint, min(y), f'x={breakpoint:.2f}', color='red', verticalalignment='bottom')
    plt.title('Original STS Fitting')
    plt.tick_params(direction='in')
    plt.legend()
    plt.show()

    x_values = x_hat
    slope_1,slope_2, slope_3 = my_pwlf.slopes[0], my_pwlf.slopes[1], my_pwlf.slopes[2]
    intercept_1, intercept_2 = my_pwlf.intercepts[0], my_pwlf.intercepts[2]

    junction_y_value_2 = my_pwlf.predict([breaks[2]])[0]
    junction_y_value_1 = my_pwlf.predict([breaks[1]])[0]
    new_break_point = ((junction_y_value_2+junction_y_value_1)/2-intercept_1)/slope_1

    y_values = np.piecewise(x_values, [x_values <= new_break_point, (x_values > new_break_point) & (x_values < breaks[2]), x_values >= breaks[2]],
                            [lambda x: slope_1 * x + intercept_1,   # 第一段直線
                            lambda x: (junction_y_value_2+ new_break_point*slope_1+intercept_1)/2,    # 第二段直線
                            lambda x: slope_3 * x + intercept_2])  # 第三段直線
    plt.scatter(x, y, label='Data')
    plt.plot(x_values, y_values, color='black', label='Fixed Fitting line')
    plt.axvline(x=breaks[2], color='red', linestyle='--')
    plt.text(breaks[2], min(y), f'x={breaks[2]:.2f}', color='red', verticalalignment='bottom')
    plt.axvline(x=new_break_point, color='red', linestyle='--')
    plt.text(new_break_point, min(y), f'x={new_break_point:.2f}', color='red', verticalalignment='bottom')
    plt.legend(loc = 'best')
    plt.tick_params(direction='in')
    plt.title('Fixed STS Fitting')
    plt.show()
    bandgap = f"{np.abs(new_break_point - breaks[2]):.2f}"
    return bandgap