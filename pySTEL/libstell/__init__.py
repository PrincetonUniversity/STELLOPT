def get_s123():
    import os
    import ctypes as ct
    from subprocess import Popen, PIPE

    # Load Libraries
    try:
        path = os.environ["STELLOPT_PATH"]+"/LIBSTELL/Release/libstell.so"
        libstell = ct.cdll.LoadLibrary(path)
        qtCreatorPath=os.environ["STELLOPT_PATH"]
    except KeyError:
        print("Please set environment variable STELLOPT_PATH")
        sys.exit(1)
    # get the full function list
    out = Popen(args="nm "+path, 
                shell=True, 
                stdout=PIPE).communicate()[0].decode("utf-8")
    attrs = [ i.split(" ")[-1].replace("\r", "") 
              for i in out.split("\n") if " T " in i]
    func = '_write_wout_file'
    module = 'read_wout_mod_' 
    names = [ s for s in attrs if module in s and func in s]
    name = names[0].replace(module, ',')
    name = name.replace(func, ',')
    s1, s2, s3 = name.split(',')
    # This is a catch for weird Macports OSX behavior
    if s1=='___':
        s1='__'
    return s1, s2, s3
s1, s2, s3 = get_s123()

__all__ = ["libstell","stellopt","fieldlines","beams3d"]
