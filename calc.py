def lvsolve(slackV,L1P,L1Q,L2P,L2Q,L3P,L3Q,phaseconfig):
    import numpy as np
    import VUFer
    reload(VUFer)
    # this is a solution based on Berg et al and the nodal method
    # it accepts a phaseconfig array describing the phase configuration at each branch point along the main feeder.
    # e.g. [3,0,0] means that all three loads are connected to
    # It has been validated as having close agreement with PSCAD within 0.0181% of load voltage magnitudes (e.g. with Fixed apparent power loads, real power values sampled from a gamma distribution* at 0.9 power factor. Unbalanced: Selected randomly with a bias towards phase 1 (40%,30%,30%). Totals on each phase were 25, 16 and 19. ).
    nonodes=len(L1P)
    vsource= np.array([slackV,slackV *(-0.5-0.8660254037844387j),slackV *(-0.5+0.8660254037844387j)]).reshape(3,)
    Iload=np.zeros((1,3),dtype=complex)
    Vload=np.zeros((1,3),dtype=complex)

    # Initialise working array
    nodedata = np.zeros((3,13,nonodes),dtype=complex)
    S= np.ones((3,nonodes),dtype=complex)
    S[0,:] = L1P + L1Q*1j
    S[1,:] = L2P + L2Q*1j
    S[2,:] = L3P + L3Q*1j
    S[S==0]=0.0000001  # set any zeros to very small number to avoid divide by 0 warning and subsequent errors
    for node in range(0,nonodes):
        nodedata[:,11,node] = S[:,node] #
        if np.any(phaseconfig[node,:]==2):
            z=np.asarray(np.where(phaseconfig[node,:]==0))
            y=np.asarray(np.where(phaseconfig[node,:]==2))
            S[y,node]=S[y,node]+S[z,node]
            S[z,node] =0.0000001
        if np.any(phaseconfig[node,:]==3):
            y=np.asarray(np.where(phaseconfig[node,:]==3))
            S[y,node]=np.sum(S[:,node])
        for r in [0,1,2]:
            if r == y:
                pass
            else:
                S[r,node] =0.0000001

    # # impedances: for 32 nodes
    # zlinea=0.00154+0.000694j # line impedance of each phase
    # zlineb=0.003+0.000703j # line impedance of each phase
    # zservice = 0.0255+0.00123j
    # zna=0.00154+0.000694j# neutral impedance
    # znb=0.003+0.00015j# neutral impedance


    #impedances for 20 nodes
    zlinea=0.00246+0.00111j# line impedance of each phase
    zlineb=0.0048+0.001125j# line impedance of each phase
    zservice = 0.0255+0.00123j#


    zna=0.00246+0.00111j# neutral impedance
    znb=0.0048+0.00024j# neutral impedance


    # admittances:
    ylinea = 1/zlinea
    ylineb = 1/zlineb
    yservice = 1/zservice
    yna = 1/zna
    ynb = 1/znb
    # put initial values working array
    for node in range(0,nonodes):
        nodedata[:,7,node]=vsource
            if node >(nonodes/2)-1:
                nodedata[0,1,node] = ylineb
                nodedata[0,2,node] = ynb
            else:
                nodedata[0,1,node] = ylinea
                nodedata[0,2,node] = yna

    nodedata[1,1,node] = yservice
    nodedata[1,2,node] = yservice
    if np.any(phaseconfig[node,:]==2):
        nodedata[:,3,node] =  1/(nodedata[:,7,node]*np.conj(nodedata[:,7,node]/nodedata[:,11,node])+[2*zservice])
        nodedata[:,0,node] = 1/(nodedata[:,7,node]*np.conj(nodedata[:,7,node]/S[:,node])+[2*zservice])
        z=np.asarray(np.where(phaseconfig[node,:]==0))
        w=np.asarray(np.where(phaseconfig[node,:]==2))
        x=np.asarray(np.where(phaseconfig[node,:]==1))
        Y1 = 1/((nodedata[w,7,node]*np.conj(nodedata[w,7,node]/nodedata[w,11,node])) + [2*zservice])
        Y2 = 1/((nodedata[z,7,node]*np.conj(nodedata[z,7,node]/nodedata[z,11,node])) + [2*zservice])
        nodedata[w,0,node] = Y1 + Y2
        nodedata[z,0,node] = 0
        nodedata[w,12,node] = 1/(nodedata[w,7,node]*np.conj(nodedata[w,7,node]/nodedata[w,11,node]))
        nodedata[z,12,node] = 1/(nodedata[w,7,node]*np.conj(nodedata[w,7,node]/nodedata[z,11,node]))
        nodedata[x,12,node] = 1/(nodedata[x,7,node]*np.conj(nodedata[x,7,node]/nodedata[x,11,node]))
    elif np.any(phaseconfig[node,:]==3):
        w=np.asarray(np.where(phaseconfig[node,:]==3))
        z=np.asarray(np.where(phaseconfig[node,:]!=3))
        nodedata[:,3,node] =  1/(nodedata[w,7,node]*np.conj(nodedata[w,7,node]/nodedata[:,11,node])+[2*zservice])
        nodedata[w,0,node] = 1/(nodedata[w,7,node]*np.conj(nodedata[w,7,node]/S[w,node])+[2*zservice])
        nodedata[z,0,node] = 0
        nodedata[:,12,node] = 1/(nodedata[w,7,node]*np.conj(nodedata[w,7,node]/nodedata[:,11,node]))
    else:
        nodedata[:,3,node] =  1/(nodedata[:,7,node]*np.conj(nodedata[:,7,node]/nodedata[:,11,node])+[2*zservice])
        nodedata[:,0,node] = 1/(nodedata[:,7,node]*np.conj(nodedata[:,7,node]/nodedata[:,11,node])+[2*zservice])
        nodedata[:,12,node] = 1/(nodedata[:,7,node]*np.conj(nodedata[:,7,node]/nodedata[:,11,node]))
    nodedata[:,4,node] = S[:,node]
    
    
    
    
    
    # Yloadcalc = The Ys used in the calculations, includes 'dummy' very high values to represent open circuits for certain phase configurations
    # Yloadserv the 3 Ys associated with the input S values combined with the service cable impedance
    #   Yload - the 3 Ys associated with the input S values alone
    
    
    #Program structure
    # 1. set up dpns ('driving point network' - see Berg-Hawkins-Pleines Paper
    # 2. calc vs dpns
    # 3. calc load vs
    # 4 recalc yload using calculated V and given S
    # repeat

    Vabsprev = np.ones((3,1,nonodes),dtype=complex)
    n=0
    for r in range(1,150):#150
        for k in range(0,nonodes):       # 1. set up dpns
            node = nonodes - k - 1 # from nonodes to 1
            if node != nonodes-1:
                # combine line Ys to to - node dpn - star to mesh transformations
                # remove node1
                A= nodedata[0,1,node+1] # Line Admittance
                B= nodedata[0,5,node+1] # Ydpn1
                C= nodedata[0,6,node+1] # YDdpn12
                D = nodedata[2,6,node+1] # YDdpn31
                E = nodedata[1,6,node+1] # YDdpn23
                F = nodedata[1,5,node+1] # Ydpn2
                G = nodedata[2,5,node+1] # Ydpn3

                D1 = A+B+C+D
                nodedata[2,6,node] = A*D/D1 # =Yli*YDdpn31/sum
                nodedata[0,6,node] = A*C/D1 # =Yli*YDdpn12/sum
                nodedata[0,5,node] = A*B/D1 # =Yli*Ydpn1/sum
                nodedata[1,6,node] = C*D/D1 + E # =YDdpn12*YDdpn31/sum + YDdpn23
                nodedata[1,5,node] = B*C/D1 + F # =Ydpn1*YDdpn12/sum + Ydpn2
                nodedata[2,5,node] = B*D/D1 + G # =Ydpn1*YDdpn31/sum + Ydpn3
                # remove node1
                B= nodedata[0,5,node]
                C= nodedata[0,6,node]
                D = nodedata[2,6,node]
                E = nodedata[1,6,node]
                F = nodedata[1,5,node]
                G = nodedata[2,5,node]
                # remove node2
                D2 = C+F+E+A
                nodedata[1,6,node] = A*E/D2
                nodedata[0,6,node] = C*A/D2
                nodedata[1,5,node] = A*F/D2
                nodedata[0,5,node] = F*C/D2 + B
                nodedata[2,6,node] = C*E/D2 + D
                nodedata[2,5,node] = F*E/D2 + G
                B= nodedata[0,5,node]   #Ydpn1
                C= nodedata[0,6,node]   #YDdpn12
                D = nodedata[2,6,node]   #YDdpn31
                E = nodedata[1,6,node]   #YDdpn23
                F = nodedata[1,5,node]   #Ydpn2
                G = nodedata[2,5,node]   #Ydpn3
                # remove node3
                D3 = D+E+G+A
                nodedata[2,6,node] = A*D/D3
                nodedata[1,6,node] = A*E/D3
                nodedata[2,5,node] = A*G/D3
                nodedata[0,5,node] = G*D/D3 + B
                nodedata[1,5,node] = G*E/D3 + F
                nodedata[0,6,node] = D*E/D3 + C
                # remove nodeN
                B= nodedata[0,5,node]
                C= nodedata[0,6,node]
                D = nodedata[2,6,node]
                E = nodedata[1,6,node]
                F = nodedata[1,5,node]
                G = nodedata[2,5,node]
                H = nodedata[0,2,node+1]#
                H = nodedata[0,2,node+1]
                D4 = B+F+G+H
                nodedata[0,5,node] = H*B/D4
                nodedata[1,5,node] = H*F/D4
                nodedata[2,5,node] = H*G/D4
                nodedata[0,6,node] = B*F/D4 + C
                nodedata[1,6,node] = F*G/D4 + E
                nodedata[2,6,node] = B*G/D4 + D
                # add dpn addmittances to Y addmittances of current node (as in parallel)
                nodedata[:,5,node] = nodedata[:,0,node] + nodedata[:,5,node]
                #print 'node:'+repr(node)
        else:
            nodedata[:,6,node] = 0  # No delta admittances at end node
            nodedata[:,5,node] = nodedata[:,0,node] #Ydpn admittances = Y admittances of end node



    for node in range(0,nonodes):
        # 2. calc v dpns - Nodal Analysis
        if node < nonodes-1:
            Y = np.array([[nodedata[0,5,node]+nodedata[0,1,node],-nodedata[0,6,node],-nodedata[2,6,node],-nodedata[0,5,node]],
            [-nodedata[0,6,node],nodedata[1,5,node]+nodedata[0,1,node],-nodedata[1,6,node],-nodedata[1,5,node]],
            [-nodedata[2,6,node],-nodedata[1,6,node],nodedata[2,5,node]+nodedata[0,1,node],-nodedata[2,5,node]],
            [-nodedata[0,5,node],-nodedata[1,5,node],-nodedata[2,5,node],nodedata[0,2,node]+np.sum(nodedata[:,5,node])]])

            Iin = np.append(nodedata[:,7,node],nodedata[2,10,node])*np.append([nodedata[0,1,node]]*3,nodedata[0,2,node]) #Iin = Vin*Yin [Zline Zneut]
            V = np.dot(np.linalg.inv(Y),Iin) ## reverted back (from above line) 28/3/2014
            nodedata[:,7,node+1]= V[0:3] ###
            nodedata[2,10,node+1]= V[3]
        else:
            Y = np.array([[nodedata[0,5,node]+nodedata[0,1,node],-nodedata[0,6,node],-nodedata[2,6,node],-nodedata[0,5,node]],
            [-nodedata[0,6,node],nodedata[1,5,node]+nodedata[0,1,node],-nodedata[1,6,node],-nodedata[1,5,node]],
            [-nodedata[2,6,node],-nodedata[1,6,node],nodedata[2,5,node]+nodedata[0,1,node],-nodedata[2,5,node]],
            [-nodedata[0,5,node],-nodedata[1,5,node],-nodedata[2,5,node],nodedata[0,2,node]+np.sum(nodedata[:,5,node])]])

            Iin = np.append(nodedata[:,7,node],nodedata[2,10,node])*np.append([nodedata[0,1,node]]*3,nodedata[0,2,node])
            V = np.dot(np.linalg.inv(Y),Iin)

      # get phase currents
        Vtrunkn = V[0:3]-V[3]
        Iphase = Vtrunkn*nodedata[:,0,node]

        #get load currents (current dividers where unbalanced phase configuration)
        if np.any(phaseconfig[node,:]==2):
            z=np.asarray(np.where(phaseconfig[node,:]==0))
            w=np.asarray(np.where(phaseconfig[node,:]==2))
            x=np.asarray(np.where(phaseconfig[node,:]==1))
            Vload = Vload.reshape((1,3))
            Iload[:,x]  = Iphase[x]
            Iload[:,w]= nodedata[w,3,node]*(Iphase[w])/(nodedata[w,3,node]+nodedata[z,3,node])
            Iload[:,z]  = nodedata[z,3,node]*(Iphase[w])/(nodedata[w,3,node]+nodedata[z,3,node])
            Vload[:,x] = (1/nodedata[x,12,node])*Vtrunkn[x]/((1/nodedata[x,12,node])+2*zservice)
            Vload[:,w] = (1/nodedata[w,12,node])*Vtrunkn[w]/((1/nodedata[w,12,node])+2*zservice)
            Vload[:,z] = (1/nodedata[z,12,node])*Vtrunkn[w]/((1/nodedata[z,12,node])+2*zservice)
        elif np.any(phaseconfig[node,:]==3):
            w=np.asarray(np.where(phaseconfig[node,:]==3))
            Iload  = nodedata[:,3,node]*Iphase[w]/(np.sum(nodedata[:,3,node]))
            Vload  = nodedata[:,11,node]/np.conj(Iload[:])
            Vload  = Iload[:]/nodedata[:,12,node]
            Vload = (1/nodedata[:,12,node])*Vtrunkn[w]/((1/nodedata[:,12,node])+2*zservice)
        else:
            Iload  = Iphase.reshape((1,3))
            Vload  = Iload[:]/nodedata[:,12,node]
            Vload = (1/nodedata[:,12,node])*Vtrunkn/((1/nodedata[:,12,node])+2*zservice)
        nodedata[:,8,node] = Iload ###
        nodedata[:,9,node] = Vload ###
        nodedata[0,10,node]= V[3]
        if node < nonodes-1:
            nodedata[1,10,node]=np.sum(Iphase) + nodedata[1,10,node+1]#np.sum(nodedata[:,8,node]) + nodedata[1,10,node+1]
        else:
            nodedata[1,10,node]=np.sum(Iphase)#np.sum(nodedata[:,8,node])

    # Check for convergence of voltage values:
    if np.all(np.abs((np.absolute(nodedata[:,9,:])-np.absolute(Vabsprev[:,0,:]))) <0.0001):#0.0001
        for node in range(0,nonodes):
            nodedata[2,1,node] = VUFer.VUFer(nodedata[0,7,node],nodedata[1,7,node],nodedata[2,7,node])
        print 'converged in ' + str(r+1) + ' iterations'
        #   print nodedata.dtype
        return nodedata
    else:
        Vabsprev[:,0,:] = nodedata[:,9,:]
        n=n+1
        ###print 'iteration' + str(r)
        if r > 148:
        ##print "did not converge"
        return 'ERROR'

    # 3 recalc yload (only gets here if no convergence)
    for node in range(0,nonodes):
        if np.any(phaseconfig[node,:]==2):
            z=np.asarray(np.where(phaseconfig[node,:]==0))
            w=np.asarray(np.where(phaseconfig[node,:]==2))
            x=np.asarray(np.where(phaseconfig[node,:]==1))
            nodedata[:,12,node] = 1/(nodedata[:,9,node]*np.conj(nodedata[:,9,node]/nodedata[:,11,node]))
            Vtrunkn = nodedata[:,7,node] - nodedata[0,10,node]
            nodedata[:,3,node] = ((nodedata[:,12,node])/(1 + 2*nodedata[:,12,node]*zservice))
            Y1 = 1/((nodedata[w,9,node]*np.conj(nodedata[w,9,node]/nodedata[w,11,node])) + [2*zservice])
            Y2 = 1/((nodedata[z,9,node]*np.conj(nodedata[z,9,node]/nodedata[z,11,node])) + [2*zservice])
            nodedata[w,0,node] = nodedata[w,3,node] + nodedata[z,3,node]
            nodedata[z,0,node] = 0
            nodedata[x,0,node] = nodedata[x,3,node]
        elif np.any(phaseconfig[node,:]==3):
            w=np.asarray(np.where(phaseconfig[node,:]==3))
            z=np.asarray(np.where(phaseconfig[node,:]!=3))
            Vtrunkn = nodedata[:,7,node]-nodedata[0,10,node]
            nodedata[:,12,node] = 1/(nodedata[:,9,node]*np.conj(nodedata[:,9,node]/nodedata[:,11,node]))
            nodedata[:,3,node] = ((nodedata[:,12,node])/(1 + 2*nodedata[:,12,node]*zservice))
            nodedata[z,0,node] = 0
            nodedata[w,0,node] = np.sum(nodedata[:,3,node])
        else:
            nodedata[:,12,node] = 1/(nodedata[:,9,node]*np.conj(nodedata[:,9,node]/nodedata[:,11,node]))
            nodedata[:,3,node] = ((nodedata[:,12,node])/(1 + 2*nodedata[:,12,node]*zservice))
            nodedata[:,0,node] = nodedata[:,3,node]

# WORKING ARRAY STRUCTURE (nodedata):
#       0           1       2       3          4       5       6         7         8       9       10       11          12
#0 [    Yloadcalc1  Yli     Yli-n   Yloadserv1 Sload1  Zdpn1   ZDdpn12   vsource1  iload1  Vload1  Vneut    Sloadact1   Yload1
#1      Yloadcalc2  Yserv   Yserv-n Yloadserv2 Sload2  Zdpn2   ZDdpn23   vsource2  iload2  Vload2  Ineut    Sloadact2   Yload2
#2      Yloadcalc3  0       0       Yloadserv3 Sload3  Zdpn3   ZDdpn31   vsource3  iload3  Vload3  VsourceN Sloadact3   Yload3]
