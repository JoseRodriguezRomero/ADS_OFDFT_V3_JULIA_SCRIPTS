import sys

if float(sys.argv[3]) % 2 == 0:
    file_names = [
        "GeometryOptimizations/singlet_" + sys.argv[1] + sys.argv[2] + ".in",
        "GeometryOptimizations/triplet_" + sys.argv[1] + sys.argv[2] + ".in",
        "GeometryOptimizations/quintet_" + sys.argv[1] + sys.argv[2] + ".in"
    ]
    s0 = 1
else:
    file_names = [
        "GeometryOptimizations/doublet_" + sys.argv[1] + sys.argv[2] + ".in",
        "GeometryOptimizations/quartet_" + sys.argv[1] + sys.argv[2] + ".in",
        "GeometryOptimizations/sextet_" + sys.argv[1] + sys.argv[2] + ".in"
    ]
    s0 = 0

for i in range(len(file_names)):
    fileID = open(file_names[i], "w")
    fileID.write("molecule { \n")
    fileID.write("-1 " + str(2*(i + 1) - s0) + "\n")
    fileID.write(sys.argv[1] + "   0.00000        0.60000        0.00000\n")
    fileID.write(sys.argv[2] + "   0.00000       -0.60000        0.00000\n")
    fileID.write("}\n\n")
    
    fileID.write("set { \n")
    fileID.write("basis              cc-pvdz\n")
    fileID.write("reference          uks\n")
    fileID.write("maxiter            500\n")
    fileID.write("geom_maxiter       500\n")
    fileID.write("DAMPING_PERCENTAGE 1\n")
    fileID.write("guess              SAP\n")
    fileID.write("}\n\n")
    
    fileID.write("optimize('b3lyp')\n")
    fileID.close()
