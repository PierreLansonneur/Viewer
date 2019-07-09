def Clear_axes(ClearAll = False):

	if (w4.get()>0)or((ROI_open==True)and(ROI_show==True)) or ((dosi_open == True)and(isodose_show == True)):
		if(w1_moved == True)or(w4.get()>0):
			for ax in [ax1, ax4]: 
				ax.clear()
				ax.set_xticks([])
				ax.set_yticks([])

                        try:    cbaxes.clear()
                        except Exception:       pass

		if(w2_moved == True)or(w4.get()>0):
			for ax in [ax2, ax5]: 
				ax.clear()
				ax.set_xticks([])
				ax.set_yticks([])
		if(w3_moved == True)or(w4.get()>0):
			for ax in [ax3, ax6]: 
				ax.clear()
				ax.set_xticks([])
				ax.set_yticks([])
	if (ClearAll == True):
		for ax in [ax1, ax2, ax3, ax4, ax5, ax6]: 
			ax.clear()
			ax.set_xticks([])
			ax.set_yticks([])
                try:    cbaxes.clear()
                except Exception:       pass

def Update_all():
	"""
	Update every axis
	"""	
	global w1_moved, w2_moved, w3_moved
	w1_moved = w2_moved = w3_moved = True
	Update()
	w1_moved = w2_moved = w3_moved = False

def Update(): 
	"""
	Update the informations on each axes
	"""
	global im1,im2,im3,dos1,dos2,dos3
	start = time.time()
	Get_axes_lim()
	Clear_axes((w1_moved)and(w2_moved)and(w3_moved))

	if(len(np.shape(volume))==3): ### 3D files #############

		### CT
		slice_ax1 = w1.get()
		slice_ax2 = w2.get()
		slice_ax3 = w3.get()

		im1 = volume[slice_ax1,:,:].T
		im2 = volume[:,slice_ax2,:].T
		im3 = volume[:,:,slice_ax3].T

		ext7=[origin[1],origin[1]-dim_y*spacing[1],origin[2]-dim_z*spacing[2],origin[2]]
		ext8=[origin[0],origin[0]+dim_x*spacing[0],origin[2]-dim_z*spacing[2],origin[2]]
		ext9=[origin[0],origin[0]+dim_x*spacing[0],origin[1]-dim_y*spacing[1],origin[1]]

		if rot1:	im1, ext7 = np.rot90(im1), [ext7[3],ext7[2],ext7[0],ext7[1]]
		if rot2:	im2, ext8 = np.rot90(im2), [ext8[3],ext8[2],ext8[0],ext8[1]]
		if rot3:	im3, ext9 = np.rot90(im3), [ext9[3],ext9[2],ext9[0],ext9[1]]

		if w1_moved or w4.get()>0:	ax1.imshow(im1, extent=ext7)
		if w2_moved or w4.get()>0:	ax2.imshow(im2, extent=ext8)
		if w3_moved or w4.get()>0:	ax3.imshow(im3, extent=ext9)


		'''### demo
		if profile_show and ID=='102771':
			#profile_ax.plot([xp1, xp2],[yp1, yp2], 'bo-')
                        ax1.plot([-94, -3],[50, -14], 'bo-')
                        ax2.plot([65, 65],[0, -26], 'bo-') 
                        ax3.plot([-46, 12],[75, 75], 'bo-')
                        ax3.plot([-100, 9],[-2, -27], 'bo-')

		if profile_show and ID=='12420':
			#profile_ax.plot([xp1, xp2],[yp1, yp2], 'bo-')
                        ax1.plot([-29, 73],[23, -58], 'bo-')
                        ax1.plot([-62, -34],[22, 50], 'bo-') 
                        ax1.plot([-55, 1],[-56, -72], 'bo-')
                        ax3.plot([39, 5],[38, 41], 'bo-')
                        #ax3.plot([-8, 48],[42, 38], 'bo-')
                '''
		slice_ax4 = int((np.shape(im1)[1]-1)*(1 - 0.001*w4.get()))
		slice_ax5 = int((np.shape(im2)[1]-1)*(1 - 0.001*w4.get()))
		slice_ax6 = int((np.shape(im3)[1]-1)*(1 - 0.001*w4.get()))

		if w4.get()>0:
			ax1.axvline(x=ext7[1] + (ext7[0]-ext7[1])*0.001*w4.get(), ls = '-', visible = w4.get())		
			ax2.axvline(x=ext8[1] + (ext8[0]-ext8[1])*0.001*w4.get(), ls = '-', visible = w4.get())
			ax3.axvline(x=ext9[1] + (ext9[0]-ext9[1])*0.001*w4.get(), ls = '-', visible = w4.get())

			if (isodose_show == False):
				ax4.plot(im1[:,slice_ax4],np.linspace(ext7[3],ext7[2],np.shape(im1)[0]),visible=w4.get())
				ax5.plot(im2[:,slice_ax5],np.linspace(ext8[3],ext8[2],np.shape(im2)[0]),visible=w4.get())
				ax6.plot(im3[:,slice_ax6],np.linspace(ext9[3],ext9[2],np.shape(im3)[0]),visible=w4.get())

		### RP
		if RP_crop_show:
                        try:
                                UpdateRPBoundingBox(ext7, ext8, ext9, rot1, rot2, rot3, ax1)
                                UpdateRPBoundingBox(ext7, ext8, ext9, rot1, rot2, rot3, ax2)
                                UpdateRPBoundingBox(ext7, ext8, ext9, rot1, rot2, rot3, ax3)

		        except Exception:
                                        try:
                                        	del ax1.patches[-1]
                                        	del ax2.patches[-1]
                                        	del ax3.patches[-1]
                                        except Exception:	pass

		### Dosimetry
		if dosi_open and isodose_show:	

			dos1 = dosi[slice_ax1_dosi(slice_ax1)-1,:,:].T
			dos2 = dosi[:,slice_ax2_dosi(slice_ax2)-1,:].T
			dos3 = dosi[:,:,slice_ax3_dosi(slice_ax3)-1].T

			cont1_dosi=[origin_dosi[1],origin_dosi[1]-dim_y_dosi*spacing_dosi[1],origin_dosi[2],origin_dosi[2]-dim_z_dosi*spacing_dosi[2]]
			cont2_dosi=[origin_dosi[0],origin_dosi[0]+dim_x_dosi*spacing_dosi[0],origin_dosi[2],origin_dosi[2]-dim_z_dosi*spacing_dosi[2]]
			cont3_dosi=[origin_dosi[0],origin_dosi[0]+dim_x_dosi*spacing_dosi[0],origin_dosi[1],origin_dosi[1]-dim_y_dosi*spacing_dosi[1]]

			dos4 = dos1[:,slice_ax2_dosi(slice_ax4)-1]
			dos5 = dos2[:,slice_ax1_dosi(slice_ax5)-1]
			dos6 = dos3[:,slice_ax1_dosi(slice_ax6)-1]

			if rot1:	dos1,cont1_dosi,dos4 = np.rot90(dos1), [cont1_dosi[2],cont1_dosi[3],cont1_dosi[1],cont1_dosi[0]], np.rot90(dos1)[:,slice_ax3_dosi(slice_ax4)-1]
			if rot2:	dos2,cont2_dosi,dos5 = np.rot90(dos2), [cont2_dosi[2],cont2_dosi[3],cont2_dosi[1],cont2_dosi[0]], np.rot90(dos2)[:,slice_ax3_dosi(slice_ax5)-1]
			if rot3:	dos3,cont3_dosi,dos6 = np.rot90(dos3), [cont3_dosi[2],cont3_dosi[3],cont3_dosi[1],cont3_dosi[0]], np.rot90(dos3)[:,slice_ax2_dosi(slice_ax6)-1]

                        if(Dvar.get()==1):
			        if w1_moved or w4.get()>0:	cs1=ax1.contour(dos1, np.nanmax(dos1)*levels, cmap = dosemap, linewidths=1, extent=cont1_dosi)
			        if w2_moved or w4.get()>0:	cs2=ax2.contour(dos2, np.nanmax(dos2)*levels, cmap = dosemap, linewidths=1, extent=cont2_dosi)
			        if w3_moved or w4.get()>0:	cs3=ax3.contour(dos3, np.nanmax(dos3)*levels, cmap = dosemap, linewidths=1, extent=cont3_dosi)

                        if(Dvar.get()==2):
			        if w1_moved or w4.get()>0:	cs1=ax1.contour(dos1, coeff_D_PTV*levels, cmap = dosemap, linewidths=1, extent=cont1_dosi)
			        if w2_moved or w4.get()>0:	cs2=ax2.contour(dos2, coeff_D_PTV*levels, cmap = dosemap, linewidths=1, extent=cont2_dosi)
			        if w3_moved or w4.get()>0:	cs3=ax3.contour(dos3, coeff_D_PTV*levels, cmap = dosemap, linewidths=1, extent=cont3_dosi)

                        if(Dvar.get()==1 or 2):
			        try:
				        #for g in range(len(cs1.collections)):	cs1.collections[g].set_label('{0} Gy'.format(int(100*levels[g])))
				        for g in range(len(cs1.collections)):	cs1.collections[g].set_label('{0} %'.format(int(100*levels[g])))
				        if(w1_moved and Dvar.get()==1):	leg = ax1.legend(frameon=False,title="$D_{max}$")
				        if(w1_moved and Dvar.get()==2):	leg = ax1.legend(frameon=False,title="$D_{PTV}$")
				        if inv_scale and w1_moved:	
					        for text in leg.get_texts():	text.set_color("white")
                                                leg.get_title().set_color("white")
			        except Exception:	pass

                        if(Dvar.get()==3):
		                if w1_moved or w4.get()>0:      cs1=ax1.imshow(np.flipud(np.ma.masked_where(dos1<0.05*np.nanmax(dosi),dos1)), alpha=0.5, cmap = dosemap, extent=cont1_dosi)
		                if w2_moved or w4.get()>0:      cs2=ax2.imshow(np.flipud(np.ma.masked_where(dos2<0.05*np.nanmax(dosi),dos2)), alpha=0.5, cmap = dosemap, extent=cont2_dosi)
		                if w3_moved or w4.get()>0:      cs3=ax3.imshow(np.flipud(np.ma.masked_where(dos3<0.05*np.nanmax(dosi),dos3)), alpha=0.5, cmap = dosemap, extent=cont3_dosi)
                                try:    
                                        cbar = P.colorbar(cs1, cax=cbaxes)
                                        if inv_scale and w1_moved:      cbaxes.tick_params(colors = "white")#,grid_color = "white")
			        except Exception:	pass

			### profiles
			ax4.plot(dos4,np.linspace(cont1_dosi[2],cont1_dosi[3],np.shape(dos1)[0]),visible=w4.get())
			ax5.plot(dos5,np.linspace(cont2_dosi[2],cont2_dosi[3],np.shape(dos2)[0]),visible=w4.get())
			ax6.plot(dos6,np.linspace(cont3_dosi[2],cont3_dosi[3],np.shape(dos3)[0]),visible=w4.get())

		### profile tool
		if profile_show:
			profile_ax.plot([xp1, xp2],[yp1, yp2], 'bo-')
			Update_profile()

		### ROI
		if (ROI_open==True)and(ROI_show==True):	
			for ROI_index in range(N_ROI):	UpdateROI(ROI_infos[ROI_index,:], slice_ax1, slice_ax2, slice_ax3, ext8, ext9) 

	"""
	if(len(np.shape(volume))==2): ### 2D files #############

		ax1.imshow(volume, extent=[origin[1],origin[1]-dim_y*spacing[1],origin[0]-dim_x*spacing[0],origin[0]])
		ax1.axvline(x=(dim_y-1)*0.001*w4.get(), ls = '-', visible = w4.get())

		if (isodose_show == False):	ax4.plot(volume[:,int((dim_y-1)*0.001*w4.get())],range(dim_x,0,-1), visible = w4.get())

		if (dosi_open == True) and (isodose_show == True):	
			cs1=ax1.contour(dos, np.nanmax(dos)*levels, cmap = dosemap,linewidths=1)
			for i in range(len(cs1.collections)):	cs1.collections[i].set_label('{0} %'.format(int(100*levels[i])))
			ax1.legend(frameon=False)
			ax4.plot(dos[:,int((dim_y-1)*0.001*w4.get())],range(dim_x,0,-1), visible = w4.get())
	"""

	t1.set_text("{0}".format( slice_ax1 ))
    	t2.set_text("{0}".format( slice_ax2 ))
    	t3.set_text("{0}".format( slice_ax3 ))

	Set_axes_lim()
	#graph1.draw() # not needed if toolbar.draw()

	if c_scale=='HOUNSFIELD':	SetDisplayRange(-1000,2000)
	if c_scale=='USER':		SetDisplayRange(0,30)

	toolbar.draw()


