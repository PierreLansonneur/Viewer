def Window_Gamma():
        """
        Display Gamma index map and passing rate
        """
	window_gamma = Toplevel(window)
	window_gamma.geometry("{0}x{1}+{2}+{3}".format(900, 435, int(0.53*w_px), 0))
	window_gamma.title('Gamma index')

	fig2 = P.figure(facecolor='lightgrey', figsize=(9, 4))
	graph2 = FigureCanvasTkAgg(fig2, master=window_gamma)
	canvas2 = graph2.get_tk_widget()
	canvas2.grid(row=0, column=0, columnspan=2, sticky=W+E+N+S)
	toolbar_frame2 = Frame(window_gamma)
	toolbar2 = NavigationToolbar2TkAgg( graph2, toolbar_frame2 )
	toolbar_frame2.grid(row=1, column=0, sticky=W+E+N+S)

	gs2 = gridspec.GridSpec(nrows=1, ncols=2, bottom=0.15, top=0.95, left = 0.07, right=0.99, hspace=0.2, wspace=0.1)
	ax7 = fig2.add_subplot(gs2[0])
	ax8 = fig2.add_subplot(gs2[1])

	ax7.set_xlabel('x (mm)')
	ax7.set_ylabel('y (mm)')
	ax8.set_xlabel('Gamma Index')

        sim1_3D, sim2_3D = np.zeros((200,200,10)), np.zeros((200,200,10))
        #sim1_3D, sim2_3D = np.zeros((200,200,40)), np.zeros((200,200,40))

        w1_get = w1.get()
        #for k,w1_set in zip(range(30),range(65,105)):
        for k,w1_set in zip(range(10),range(w1_get-5,w1_get+5)):
                w1.set(w1_set)
        	Update_all()
                sim1, sim2 = InterpolateDosiCT(ax1,-100,-100,100,100)
                #sim1, sim2 = InterpolateDosiCT(ax1,-150,-100,50,100)
                sim1_3D[:,:,k], sim2_3D[:,:,k] = sim1, sim2

        '''
        ### Normalize to the same max dose
        sim1_3D = (100./np.nanmax(sim1_3D))*sim1_3D
        sim2_3D = (100./np.nanmax(sim2_3D))*sim2_3D
        '''
        ### Normalize to the same avge dose
        sim1_3D = (100./np.mean(sim1_3D))*sim1_3D
        sim2_3D = (100./np.mean(sim2_3D))*sim2_3D

        ### Compute gamma index
        #gamma = gamma_matrix3D(sim1_3D, sim2_3D)                       # 2% / 2 mm
        gamma = GammaIndex(sim1_3D, sim2_3D, DTA=3, dmax=0.03)          # 3% / 3 mm
        #gamma = sim1_3D-sim2_3D                                        # dose difference

        a7 = ax7.imshow(gamma[:,:,5], cmap=P.get_cmap('coolwarm'))
        a7.set_clim(0, 2) 

	ax8.hist(np.ravel(gamma), bins='auto', histtype='step', density=True)
        ax8.set_yscale('log')

        graph2.draw()
        toolbar2.draw()
        fig2.colorbar(a7, ax=ax7, pad=-0.001)
        #fig2.savefig('./tmp/GammaIndex.pdf',  bbox_inches='tight')
        w1.set(w1_get)

def InterpolateDosiCT(ax,xp1,yp1,xp2,yp2):
        """
        return the matrices of CT and dosi displayed on axe 'ax'
        The 2D matrices are interpolated between the coordinates (xp1,yp1) and (xp2,yp2) with a spacing of 1 mm
        """
        # dosi matrix
        xp1_vox, yp1_vox = mm2voxel(xp1, yp1, ax, tag='dosi')
        xp2_vox, yp2_vox = mm2voxel(xp2, yp2, ax, tag='dosi')
        xp,yp = np.linspace(xp1_vox,xp2_vox,200), np.linspace(yp1_vox,yp2_vox,200)
        sim1 = map_coordinates(AxisImage(ax,tag='dosi'), np.meshgrid(xp, yp))

        # CT matrix
        xp1_vox, yp1_vox = mm2voxel(xp1, yp1, ax)
        xp2_vox, yp2_vox = mm2voxel(xp2, yp2, ax)
        xp,yp = np.linspace(xp1_vox,xp2_vox,200), np.linspace(yp1_vox,yp2_vox,200)
        sim2 = map_coordinates(AxisImage(ax,tag='CT'),np.meshgrid(xp, yp))

        return sim1, sim2

def GammaIndex(m1, m2, DTA=2, dmax=0.02):
	# Compute matrix of gammma indices.
	# m1: 	reference matrix
	# m2: 	tested matrix
	# DTA: 	maximum distance-to-agreement (in voxels)
	# dmax: maximum dose difference

        print 'Calculating Gamma index..'

	if m1.shape != m2.shape:	print "Error: cannot compute for matrices of different sizes."

	dmax = dmax*np.nanmax(m1)
	n = 0

        ### 2D matrix
	if(len(m1.shape)==2):
                tmp = np.zeros((m1.shape[0]-2*DTA,m1.shape[1]-2*DTA,4*DTA**2))

	        for i in range(-DTA,DTA):
		        for j in range(-DTA,DTA):
		                result1 = np.power( (m1[i+DTA:i-DTA,j+DTA:j-DTA] - m2[DTA:-DTA,DTA:-DTA])/dmax,2 )
		                result2 = (i*i + j*j)/DTA*DTA 
		                tmp[:,:,n] = np.sqrt(result1+result2)
		                n = n+1

                output = np.nanmin(tmp, axis = 2)
                GIPR = 100*len(np.where(output <= 1)[0])/(output.shape[0]*output.shape[1])

        ### 3D matrix
	if(len(m1.shape)==3):   
                tmp = np.zeros((m1.shape[0]-2*DTA,m1.shape[1]-2*DTA,m1.shape[2]-2,8*DTA**2))
        	#tmp = np.zeros((m1.shape[0]-2*DTA,m1.shape[1]-2*DTA,m1.shape[2]-2*DTA,8*DTA**3))

	        for i in range(-DTA,DTA):
		        for j in range(-DTA,DTA):
		                #for k in range(-DTA,DTA):
		                for k in range(-1,1):  
			                #result1 = np.power( (m1[i+DTA:i-DTA,j+DTA:j-DTA,k+DTA:k-DTA] - m2[DTA:-DTA,DTA:-DTA,k+DTA:k-DTA])/dmax,2 )
			                result1 = np.power( (m1[i+DTA:i-DTA,j+DTA:j-DTA,k+1:k-1] - m2[DTA:-DTA,DTA:-DTA,k+1:k-1])/dmax,2 )
			                #result2 = (i*i + j*j + k*k)/DTA**2 
			                result2 = (i*i + j*j + (1./9)*k*k)/DTA**2 
			                tmp[:,:,:,n] = np.sqrt(result1+result2)
			                n = n+1

	        output = np.nanmin(tmp, axis = 3)
                GIPR = 100.*len(np.where(output <= 1)[0])/(output.shape[0]*output.shape[1]*output.shape[2])

        print '{0:.2f} % of pixels below 1'.format(GIPR)

	return output

