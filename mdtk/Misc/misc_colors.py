
"""
# My color palettes for matplotlib
# usage: 
import matplotlib.pyplot as plt
from mdtk.Misc.misc_colors import ColorPalette_classic_01_12
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color=ColorPalette_classic_01_12)

"""
# Define color palettes
# Name format: ColorPalette_[subtitle]_[index]_[# of colors]
ColorPalette_aesthetic_01_5=["#66545e","#a39193","#aa6f73","#eea990","#f6e0b5"]
ColorPalette_classic_01_12=["#ffa510", "#188ac9", "#fbda41" , "#f74d4d" , "#2455a4" , "#41b7ac", "#012c56","#feaf0f","#0a8ecf","#ffca6c","#fb5050", "#002c53"]
# ColorPalette_classic_01_12=["#ffa510", "#0c84c6", "#ffbd66" , "#f74d4d" , "#2455a4" , "#41b7ac", "#012c56","#feaf0f","#0a8ecf","#ffca6c","#fb5050", "#002c53"]
# ColorPalette_classic_01_12=["#ffca6c", "#ffa510","#869d9d", "#0c84c6", "#e77c8e" , "#f74d4d" , "#2455a4" , "#41b7ac", "#012c56","#daa45a","#5c2223", "#002c53"]
ColorPalette_classic_02_5=["#8ECAE6","#219EBC","#023047","#FFB703","#FB8500"]
ColorPalette_light_01_12=["#63b2ee" , "#76da91" , "#f8cb7f" , "#f89588" , "#7cd6cf" , "#9192ab" , "#7898e1" , "#efa666" , "#eddd86" , "#9987ce" , "#63b2ee" , "#76da91"] 
ColorPalette_dark_01_12=["#3b6291", "#943c39", "#779043", "#624c7c", "#388498", "#bf7334", "#3f6899", "#9c403d", "#7d9847", "#675083", "#3b8ba1", "#c97937"]
ColorPalette_pale_01_10=["#E2E2DF","#D2D2CF","#E2CFC4","#F7D9C4","#FAEDCB","#C9E4DE","#C6DEF1","#DBCDF0","#F2C6DE","#F9C6C9"]

# list all color palettes
def list_all_color_palettes():
    # get all global variables
    global_vars = globals()
    # use regular expression to find all color palettes, ColorPalette_[subtitle]_[index]_[# of colors]
    import re
    pattern = re.compile(r"ColorPalette_\w+_\d+_\d+")
    ColorPalettes = [var for var in global_vars if pattern.match(var)]
    for i, cp in enumerate(ColorPalettes):
        if i < len(ColorPalettes)-1:
            print(cp, end="; ")
        else:
            print(cp)