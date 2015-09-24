# -*- coding: utf-8 -*-

import os, sys

def dll2lib(dll):

    txt = dll.replace(".dll",".txt")
    ddef = dll.replace(".dll",".def")
    lib = dll.replace("lib","").replace(".dll",".lib")

    os.system("dumpbin /nologo /exports /out:%s %s"%(txt,dll))

    txtfile = open(txt,"r")
    txtcont = txtfile.readlines()

    txtfile.close()

    deffile = open(ddef,"w")
    deffile.write("LIBRARY %s\n"%dll)
    deffile.write("EXPORTS\n")

    n = 16
    cl = txtcont[n]

    while cl!="\n":
        cl_split = cl.split()
        func_order = cl_split[0]
        if "=" in cl:
            func_name = cl_split[-3]
        else:
            func_name = cl_split[-1]
        deffile.write("\t%s @%s\n"%(func_name,func_order))
        n += 1
        cl = txtcont[n]

    deffile.close()
    os.remove(txt)

    os.system("lib /nologo /def:%s /machine:x64 /out:%s"%(ddef,lib))
    os.remove(lib.replace("lib","exp"))
    # os.remove(ddef)

if __name__ == '__main__':
    for dll in sys.argv[1:]:
        dll2lib(dll)