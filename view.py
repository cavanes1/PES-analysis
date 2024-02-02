# parameters
rnge = 2 # x-axis range in each direction for a plot
step = 0.005 # Step size, 0.04 is not bad
origin_name = "C2vCI" # Path of initial Cartesian coordinate geometry file

# this program was built based off surfcurve.py as a starting point
# cartgrd write and energy.all write need to be edited for the molecule of interest
# Before running, make sure you have
#   [origin_name]
#   intcfl (required by COLUMBUS's cart2int)
#   fit.in
#   dat.x, which itself additionally requires:
#       geom.all
#       basis.in
#       Hd.CheckPoint
# This program automatically generates the required files
#   energy.all
#   names.all
#   refgeom
# Output: plots in PLOTS named PLOTS/*.png

# module import
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import os
import subprocess
print("\nModules imported")

# read fit.in
nstates = 0
h = open("fit.in", "r")
fitin = h.readlines()
h.close()
for param in fitin:
    if "eshift" in param:
        eshift = float(param.split()[-1][:-3])
    elif "nstates" in param:
        nstates = int(param.split()[-1])
print("eshift  = " + str(eshift))
print("nstates = " + str(nstates))
conversion = 219474 # cm-1 per hartree

# generate cart2intin
c2itxt = """ &input
    calctype='cart2int'
 /"""
i2ctxt = """ &input
    calctype='int2cart'
 /"""
f = open("./cart2intin", "w")
f.write(c2itxt)
f.close()

# generate dummy cartgrd
cartgrdtxt = """    0.739462D-07   0.107719D-01   0.359672D-02
  -0.637601D-07   0.573448D-03   0.113684D-01
  -0.408939D-07  -0.123975D-01  -0.519692D-02
   0.345788D-07  -0.168148D-02  -0.133691D-01
  -0.145659D-07   0.423517D-02   0.559146D-02
   0.106949D-07  -0.150151D-02  -0.199053D-02"""
f = open("./cartgrd", "w")
f.write(cartgrdtxt)
f.close()

# run cart2int to convert to internals
os.system("cp " + origin_name + " geom")
os.system("cp " + origin_name + " origgeom")
rv = subprocess.run(["/home/cavanes1/col/Columbus/cart2int.x"],cwd="./",capture_output=True)
print("\nConverted Cartesian origin geometry to internal coordinates")
f = open("./cart2intin", "w")
f.write(i2ctxt)
f.close()

# read given internal coordinate geometry file
f = open("intgeom", "r")
init = f.readlines()
f.close()

# extract data from internal coordinate geometry
init_data = []
for line in init:
    init_data.append(float(line))
init_data = np.array(init_data)
intdim = len(init_data)
print("\nInternal coordinate geometry data read and extracted")

# generate curve endpoints
endpts = []
names = []
sr = np.sqrt(rnge)
# single-coordinate displacements
for i in range(1,intdim + 1):
    names.append("coord" + str(i))
    gp = init_data.copy()
    gp[i-1] += rnge
    endpts.append(gp)
# positive pairwise displacements
for i in range(1,intdim + 1):
    break # comment out
    for j in range(i,intdim + 1):
        names.append("i" + str(i) + "jplus" + str(j))
        gp = init_data.copy()
        gp[i-1] += sr
        gp[j-1] += sr
        endpts.append(gp)
# opposite pairwise displacements
for i in range(1,intdim + 1):
    break # comment out
    for j in range(i,intdim + 1):
        names.append("i" + str(i) + "jmin" + str(j))
        gm = init_data.copy()
        gm[i-1] -= sr
        gm[j-1] -= sr
        endpts.append(gm)
print("\nCurve endpoints generated\n")

# analyze present files
mypath = "./"
allfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

# back up files which will be appended
filelist = ["geom.all", "names.all", "energy.all"]
for fl in filelist:
    if fl + ".old" not in allfiles: # if backup does not exist
        os.system("mv " + fl + " " + fl + ".old")
    else:
        print("Note: Backup already exists for " + fl)
        if fl in allfiles:
            os.system("rm " + fl)
        else:
            print("    Note: " + fl + " not present")
print("Backup complete for files to be appended")

# run this code for each step
def slurmcop(target, interpvec):
    str_target = str(target)
    #print("Working on point #" + str_target)
    # generate interpolated geometry as intgeomch
    f = open("./intgeomch", "w")
    currgeom = init_data + interpvec*step*target
    for coord in currgeom:
        formcoord = format(float(coord), "14.8f") + "\n"
        f.write(formcoord)
    f.close()
    # run cart2int to convert to Cartesians
    rv = subprocess.run(["/home/cavanes1/col/Columbus/cart2int.x"],cwd="./",capture_output=True)
    # read newly-generated Cartesian geometry file
    f = open("./geom.new", "r")
    cart = f.readlines()
    f.close()
    # write to geom.all
    f = open("./geom.all", "a")
    for line in cart:
        f.write(line)
    f.close()
    # use adjacent reference geometry
    os.system("cp geom.new geom")
    # write to names.all
    f = open("./names.all", "a")
    f.write(format(target, "5d") + "   " + str_target + "\n")
    f.close()
    # write to energy.all
    f = open("./energy.all", "a")
    eallstring = "-256.787847331183  "*(nstates - 1)
    f.write(eallstring + "-256.787847331183\n")
    f.close()

# run this for each curve
def curve(name, endpoint):
    # create interpolation vector
    interpvec = (endpoint - init_data)*2
    print("\nInterpolation vector created for " + name + "\n")
    # run the above function automatically
    os.system("cp origgeom geom")
    print("Doing points on the right")
    for i in range(int(1/2/step) + 1): # 0 to 100, previously 100 to 200
        slurmcop(i, interpvec) # was i + 1
    os.system("cp origgeom geom")
    print("Doing points on the left")
    for i in range(-1, -int(1/2/step) - 1, -1): # -100 to -1, previously 99 to 0
        slurmcop(i, interpvec)

    # run dat.x
    rv = subprocess.run(["./dat.x"],cwd="./",capture_output=True)
    f = open("./dat.log", "w")
    f.write(rv.stdout.decode('utf8'))
    f.close()
    print("\nFinished running dat.x")

    # prepare surface data for plot
    xs = []
    unadjusted_energies = []
    for state in range(nstates):
        unadjusted_energies.append([])

    f = open("allenergies.csv", "r")
    lines = f.readlines()
    f.close()
    for line in lines:
        eners = line.split(",")
        xs.append(int(eners[0].strip(' \n\t')))
        for state in range(nstates):
            unadjusted_energies[state].append(float(eners[state + 1].strip(' \n\t')))
    xs = np.array(xs)
    adj_E = np.array(unadjusted_energies)
    adj_E = (adj_E + eshift)*conversion

    # make plot
    for i in range(nstates):
        #plt.ylim(-5000,90000)
        #plt.xlim(-1.5,1.5)
        plt.title(name)
        plt.xlabel('LST progression from ' + origin_name + ' (arb. u.)')
        plt.ylabel('Relative total electronic energy (cm$^{-1}$)')
        z = sorted(zip(xs*rnge*2*step, adj_E[i]))
        x=[i[0] for i in z]
        y=[i[1] for i in z]
        plt.plot(x, y, label = str(i + 1))
    plt.savefig("PLOTS/" + name)
    plt.close()
    print("\nPlot saved to " + name + ".png")

    # clean files for next plot
    for fl in filelist:
        os.system("rm " + fl)

os.system("mkdir PLOTS")
# for each plot
for i in range(len(names)):
    curve(names[i], endpts[i])

# restore backups
for fl in filelist:
    os.system("mv " + fl + ".old " + fl)
print("\nBackup restored for appended files")
