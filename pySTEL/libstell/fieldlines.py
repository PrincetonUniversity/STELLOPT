def read_fieldlines(file):
    import h5py
    import numpy as np
    fieldline_data = {}
    #read file
    with h5py.File(file,'r') as f:
        # Logicals
        for temp in ['ladvanced', 'laxis_i', 'lcoil', 'lmgrid', 'lmu', 'lpies', 'lreverse', 'lspec', 'lvac', 'lvessel', 'lvmec']:
            if temp in f:
                fieldline_data[temp] = np.int(f[temp][0])
        # Integers
        for temp in ['nlines', 'nphi', 'npoinc', 'nr', 'nsteps', 'nz']:
            if temp in f:
                fieldline_data[temp] = np.int(f[temp][0])
        # Floats
        for temp in ['VERSION','iota0']:
            if temp in f:
                fieldline_data[temp] = np.float(f[temp][0])
        # Arrays
        for temp in ['phiaxis', 'raxis', 'zaxis', 'B_lines', 'PHI_lines', 'R_lines', 'Z_lines', 'B_PHI', 'B_R', 'B_Z',\
                    'wall_vertex', 'wall_faces', 'wall_strikes', 'A_R', 'A_PHI', 'A_Z', 'L_lines', 'Rhc_lines', 'Zhc_lines']:
            if temp in f:
                fieldline_data[temp] = np.array(f[temp][:])
    # Make derived arrays
    fieldline_data['X_lines'] = fieldline_data['R_lines']*np.cos(fieldline_data['PHI_lines'])
    fieldline_data['Y_lines'] = fieldline_data['R_lines']*np.sin(fieldline_data['PHI_lines'])
    brtemp = np.zeros(fieldline_data['B_R'].shape)
    bztemp = np.zeros(fieldline_data['B_Z'].shape)
    for i in range(fieldline_data['nr']):
        brtemp[i,:,:] = fieldline_data['B_R'][i,:,:]*fieldline_data['B_PHI'][i,:,:]/fieldline_data['raxis'][i]
        bztemp[i,:,:] = fieldline_data['B_Z'][i,:,:]*fieldline_data['B_PHI'][i,:,:]/fieldline_data['raxis'][i]
    fieldline_data['B_R'] = brtemp
    fieldline_data['B_Z'] = bztemp
    return fieldline_data

def calc_iota(data):
    import numpy as np
    x = np.zeros(data['R_lines'].shape)
    y = np.zeros(data['Z_lines'].shape)
    nstep = data['nsteps']
    for i in range(data['nsteps']):
        x[i,:] = data['R_lines'][i,:]-data['R_lines'][i,0]
        y[i,:] = data['Z_lines'][i,:]-data['Z_lines'][i,0]
    theta  = np.arctan2(y,x)
    dtheta = np.diff(theta,axis=0)
    dtheta = np.where(dtheta<-np.pi,dtheta+2*np.pi,dtheta)
    dtheta = np.where(dtheta>np.pi,dtheta-2*np.pi,dtheta)
    theta  = np.cumsum(dtheta,axis=0)
    out         = {}
    out['rho']  = np.mean(np.sqrt(x*x+y*y),axis=0)
    out['iota'] = np.zeros(out['rho'].shape)
    out['iota_err'] = np.zeros(out['rho'].shape)
    for i in range(data['nlines']):
        p, residuals, rank, singular_values, rcond = np.polyfit(data['PHI_lines'][0:nstep-1,i],theta[:,i],1,full=True)
        out['iota'][i]=p[0]
        out['iota_err'][i] = np.sqrt(residuals)
    out['iota'][0]=2*out['iota'][1]-out['iota'][2]
    return out
