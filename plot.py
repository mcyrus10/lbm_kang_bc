import matplotlib.pyplot as plt
from numpy import meshgrid,array
#------------------------------------------------------------------------------
# this reads in the fields and respective values from the input file that the
# simulation used
#------------------------------------------------------------------------------
def readPalabos(filename,
                nx,
                ny):
    with open(filename,'r') as f:
        text = array(f.read().split(" "))[0:-1].reshape(nx,ny)
    return(text)


def read_xml(   f_name,
                fields):
    retDict = {}
    with open(f_name,'r') as f:
        text = f.read().split("\n")
    for line in text:
        for field in fields:
            if "<{}>".format(field[0]) in line:
                temp = list(filter(None,line.split(" ")))
                retDict[field[0]] = field[1](temp[1])
    return(retDict)

fields = [
            ('nx',int),
            ('ny',int),
            ('Da',float),
            ('D',float),
            ('C_0',float),
            ('C_eq',float),
        ]
inputParameters = read_xml("params.xml",fields)

Da = inputParameters['Da']
nx = inputParameters['nx']
ny = inputParameters['ny']
D = inputParameters['D']
C0 = inputParameters['C_0']
Ceq = inputParameters['C_eq']

jx = ny
ix = nx
kr = Da*D/ny

lbm = readPalabos('lbm_contour.dat',nx,ny)

#==============================================================================
# Contour plotting (saved as .png)
#==============================================================================
fig, ax = plt.subplots()
linewidths = .5
levels = [0.1*x for x in range(10)]
X,Y = meshgrid([range(0,nx)],[range(0,ny)])
CS_lbm = ax.contour(X,Y,lbm.T,
                    linewidths = linewidths,
                    colors = ('b'),
                    levels = levels)
ax.clabel(CS_lbm, inline=1, fontsize =10)

plt.minorticks_on()
plt.grid(b=True, which = 'major')
plt.grid(b=True, which = 'minor')
plt.xlabel('x')
plt.ylabel('y')
ax.set_title("n$_x$ = {}; n$_y$= {}; Da = {:.2f}; D = {:.3f}; k$_r$ = {:.3f}".format(nx,ny,Da,D,kr))
plt.savefig("lbm_solution.png",dpi = 300)
plt.show()
