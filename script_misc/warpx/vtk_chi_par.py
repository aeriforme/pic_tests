import vtk
import numpy as np
from vtk.util import numpy_support
from scipy.constants import *
import numpy as np
from mpi4py import MPI
import openpmd_api as io
import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import *
matplotlib.rcParams.update({'figure.autolayout': True})

my_dpi=300
z=np.linspace(0,70,875)
y=np.linspace(-30,30,750)
x=np.linspace(-30,30,750)
lambda_SI = 0.8 * micron # wavelength
omega_SI = 2.0 * pi * c / lambda_SI
n_crit= (m_e*epsilon_0*(2*pi*c)**2)/((e*0.8e-6)**2)


comm = MPI.COMM_WORLD
series = io.Series("../diags/Fields_Particles/openpmd_%T.bp",io.Access.read_only)
if 0 == comm.rank:
	print("The Series contains {0} iterations:".format(len(series.iterations)))
	for i in series.iterations:
		print("\t {0}".format(i))
	print("")
	i = series.iterations[0]
	print("Iteration 0 contains {0} meshes:".format(len(i.meshes)))
	for m in i.meshes:
		print("\t {0}".format(m))
	print("")
	print("Iteration 0 contains {0} particle species:".format(len(i.particles)))
	for ps in i.particles:
		print("\t {0}".format(ps))
		print("With records:")
		for r in i.particles[ps]:
			print("\t {0}".format(r))
chunk_offset = [1,1,comm.rank + 1]
chunk_extent = [750, 750, 1]
for n in [3974]:
	i = series.iterations[n]
	ux_elefoam=i.meshes["ux_ele_foam"][io.Mesh_Record_Component.SCALAR]
	uy_elefoam=i.meshes["uy_ele_foam"][io.Mesh_Record_Component.SCALAR]
	uz_elefoam=i.meshes["uz_ele_foam"][io.Mesh_Record_Component.SCALAR]
	gamma_elefoam=i.meshes["gamma_ele_foam"][io.Mesh_Record_Component.SCALAR]
	B = i.meshes["B"]
	B_x = B["x"]
	B_y = B["y"]
	B_z = B["z"]
	E = i.meshes["E"]
	E_x = E["x"]
	E_y = E["y"]
	E_z = E["z"]
	ux = ux_elefoam.load_chunk(chunk_offset, chunk_extent)
	uy = uy_elefoam.load_chunk(chunk_offset, chunk_extent)
	uz = uz_elefoam.load_chunk(chunk_offset, chunk_extent)
	gamma = gamma_elefoam.load_chunk(chunk_offset, chunk_extent)
	bx=B_x.load_chunk(chunk_offset, chunk_extent)
	by=B_y.load_chunk(chunk_offset, chunk_extent)
	bz=B_z.load_chunk(chunk_offset, chunk_extent)
	ex=E_x.load_chunk(chunk_offset, chunk_extent)
	ey=E_y.load_chunk(chunk_offset, chunk_extent)
	ez=E_z.load_chunk(chunk_offset, chunk_extent)
	if 0 == comm.rank:
		print("Queued the loading of a single chunk per MPI rank from disk, "
              "ready to execute")
	series.flush()
	if 0 == comm.rank:
		print("Chunks have been read from disk")
	for i in range(comm.size):
		if i == comm.rank:
			print("Rank {} - Read chunk contains:".format(i))
			bx *= B_x.unit_SI
			by *= B_y.unit_SI
			bz *= B_z.unit_SI
			ex *= E_x.unit_SI
			ey *= E_y.unit_SI
			ez *= E_z.unit_SI
			Es=m_e**2*c**3/hbar/e
			print('non calcolata')
			print(gamma)
			data1=gamma/Es*np.sqrt((ex+c/gamma*(uy*bz-uz*by))**2+(ey+c/gamma*(uz*bx-ux*bz))**2+(ez+c/gamma*(ux*by-uy*bx))**2-((ux*ex+uy*ey+uz*ez)/gamma)**2)
			print('calcolata')
			data1=np.asarray(data1[::2,::2,::2])
			#data2=np.asarray(data2[::2,::2,::2])
			#data3=np.asarray(data3[::2,::2,::2])
	if comm.rank == 0:
    		recvbuf = np.empty((x.size,y.size,z.size))
	comm.Gather(data1, recvbuf, root=0)
	if 0 == comm.rank:
		data1=recvbuf
		print(data1) #.reshape((x.size,y.size,z.size))
		gx=x
		gy=y
		gz=z
		dx=(gx[1]-gx[0])
		dy=(gy[1]-gy[0])
		dz=(gz[1]-gz[0])
		ox=gx[0]
		oy=gy[0]
		oz=gz[0]
		imdata = vtk.vtkImageData()
		depthArray = numpy_support.numpy_to_vtk(np.ravel(data1, order='F'),  deep=True, array_type=vtk.VTK_DOUBLE)
		imdata.SetDimensions(data1.shape)
		#    imdata.SetSpacing([1,1,1])
		#    imdata.SetOrigin([0,0,0])
		imdata.SetSpacing([dx,dy,dz])
		imdata.SetOrigin([ox,oy,oz])
		imdata.GetPointData().SetScalars(depthArray)
		writer = vtk.vtkMetaImageWriter()
		outfile="chifoam_%s.mhd" % n
		writer.SetFileName(outfile)
		writer.SetInputData(imdata)
		writer.Write()
	del series


