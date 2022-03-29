'''
Sets plot style globally. 

'''
import matplotlib
from matplotlib import rcParams

rcParams['font.size'] = 18
rcParams['font.family'] = 'serif'
rcParams['xtick.major.size'] = 8
rcParams['xtick.labelsize'] = 'large'
rcParams['xtick.direction'] = "in"
rcParams['xtick.minor.visible'] = True
rcParams['xtick.top'] = True
rcParams['ytick.major.size'] = 8
rcParams['ytick.labelsize'] = 'large'
rcParams['ytick.direction'] = "in"
rcParams['ytick.minor.visible'] = True
rcParams['ytick.right'] = True
rcParams['xtick.minor.size'] = 4
rcParams['ytick.minor.size'] = 4
rcParams['xtick.major.pad'] = 10
rcParams['ytick.major.pad'] = 10
rcParams['legend.numpoints'] = 1
rcParams['mathtext.fontset'] = 'cm'
rcParams['mathtext.rm'] = 'serif'
rcParams['axes.labelsize'] = 'xx-large'
rcParams['axes.prop_cycle'] = matplotlib.cycler(color=['9F6CE6','FF984A','538050','6FADFA','7D7D7D','black'])