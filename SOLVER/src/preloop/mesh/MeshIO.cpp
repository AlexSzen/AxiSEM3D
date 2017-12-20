#include "MeshIO.h"
#include "Parameters.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "eigenc.h"
#include "Domain.h"
#include "Element.h"
#include "Mesh.h"
#include "Point.h"
#include "XMPI.h"
#include "SlicePlot.h"
#include "Quad.h"

MeshIO::MeshIO(const Mesh *mesh, const std::string &fname): mMesh(mesh), mFileName(fname) {
	//nothing
}


void MeshIO::dumpFields(const Domain &domain, const Source &source, const Parameters &par) {
	
	// ---------- <GET PARAMETERS> ----------
	double rmin = par.getValue<double>("RMIN")*1e3;
	double rmax = par.getValue<double>("RMAX")*1e3;
	double tmin = par.getValue<double>("THETA_MIN")*degree;
	double tmax = par.getValue<double>("THETA_MAX")*degree;
	int gll_ani = par.getValue<int>("GLL_ANIMATION");
	// ---------- </GET PARAMETERS> ----------

	// ---------- <GENERAL> ----------
	
	
	int gllpoints_proc = domain.getNumPoints();
	int gllpoints_all = XMPI::sum(gllpoints_proc);

	int elems_proc = domain.getNumElements();
	int elems_all = XMPI::sum(elems_proc);

	int elems_proc_ani = 0, elemNus_proc = 0, elemNus_proc_ani = 0;

	for (int i = 0; i < elems_proc; i++) {
		int maxNu = domain.getElement(i)->getMaxNu() + 1;
		elemNus_proc += maxNu;
		if (domain.getElement(i)->needDumping(rmin,rmax,tmin,tmax)) {
			elems_proc_ani++;
			elemNus_proc_ani += maxNu;	
		}
	}
	
	
	int elems_all_ani = XMPI::sum(elems_proc_ani);
	int elemNus_all = XMPI::sum(elemNus_proc);
	int elemNus_all_ani = XMPI::sum(elemNus_proc_ani);
	// ---------- </GENERAL> ----------
		
	// ---------- <GET FIELDS> ----------	
	std::vector<int> NusAni(elems_proc,0); 
	std::vector<int> NrsAni(elems_proc,0);
	std::vector<int> Nus(elems_proc,0); 
	std::vector<int> Nrs(elems_proc,0);
	vec_RMatPP mesh_S(elems_proc, RMatPP::Zero());
	vec_RMatPP mesh_Z(elems_proc, RMatPP::Zero());
	
	RMatXX_RM zero_matAni(gll_ani,gll_ani); //animations 
	zero_matAni.setZero();	
	vec_RMatXX_RM zero_ar2_matAni(2,zero_matAni);
	
	vec2_RMatXX_RM mesh_vp_ani;
	vec2_RMatXX_RM mesh_vs_ani;
	
	vec_ar2_RMatPP mesh_vp, mesh_vs;
	
	for (int i=0; i<elems_proc; i++) {
		
		Element* elem = domain.getElement(i);
		Quad* quad = mMesh->mQuads[i];
		 
		if (elem->needDumping(rmin,rmax,tmin,tmax)) { //for now keep the flexible dumping even for wavefields for kernels
      	NusAni[i] = elem->getMaxNu()+1;
      	NrsAni[i] = elem->getMaxNr();
		Nus[i] = elem->getMaxNu()+1;
		Nrs[i] = elem->getMaxNr();
		
		vec_CMatPP vp_elem = quad->getMaterialFourier("vp", SlicePlot::PropertyRefTypes::Property3D);
		vec_CMatPP vs_elem = quad->getMaterialFourier("vs", SlicePlot::PropertyRefTypes::Property3D);

		for (int ialpha=0; ialpha<vp_elem.size(); ialpha++) {
			
			mesh_vp_ani.push_back(zero_ar2_matAni);		
			(mesh_vp_ani.back())[0] = vp_elem[ialpha].block(nPntEdge/2-gll_ani/2, nPntEdge/2-gll_ani/2, gll_ani, gll_ani).real();
			(mesh_vp_ani.back())[1] = vp_elem[ialpha].block(nPntEdge/2-gll_ani/2, nPntEdge/2-gll_ani/2, gll_ani, gll_ani).imag();

			mesh_vs_ani.push_back(zero_ar2_matAni);		
			(mesh_vs_ani.back())[0] = vs_elem[ialpha].block(nPntEdge/2-gll_ani/2, nPntEdge/2-gll_ani/2, gll_ani, gll_ani).real();
			(mesh_vs_ani.back())[1] = vs_elem[ialpha].block(nPntEdge/2-gll_ani/2, nPntEdge/2-gll_ani/2, gll_ani, gll_ani).imag();
			
			mesh_vp.push_back(zero_ar2_RMatPP);
			(mesh_vp.back())[0] = vp_elem[ialpha].real();
			(mesh_vp.back())[1] = vp_elem[ialpha].imag();
			
			mesh_vs.push_back(zero_ar2_RMatPP);
			(mesh_vs.back())[0] = vs_elem[ialpha].real();
			(mesh_vs.back())[1] = vs_elem[ialpha].imag();


		}
		
	    } else {
		NusAni[i] = 0;
		NrsAni[i] = 0;
		Nus[i] = 0; 
		Nrs[i] = 0;
		}
		
		// coords 
		mesh_S[i] = (elem->getCoordsPoints())[0];
		mesh_Z[i] = (elem->getCoordsPoints())[1];
		
		
    }

	// ---------- </GET FIELDS> ----------	


	// ---------- <DEFINE WRITING OFFSETS> ----------
	
	int temp_startElem, temp_startElemNu;
	std::vector<int> temp_countElem(XMPI::nproc(),0);
	std::vector<int> temp_countElemNu(XMPI::nproc(),0);

	XMPI::gather(elems_proc, temp_countElem, true); // non fourier fields 
	temp_startElem = std::accumulate(temp_countElem.begin(), temp_countElem.begin()+XMPI::rank(), 0);
	
	XMPI::gather(elemNus_proc, temp_countElemNu, true); //fourier fields 
	temp_startElemNu = std::accumulate(temp_countElemNu.begin(), temp_countElemNu.begin()+XMPI::rank(),0);

	std::vector<size_t> startElem, countElem;
	startElem.push_back(temp_startElem);
	countElem.push_back(temp_countElem[XMPI::rank()]);
	
	std::vector<size_t> startElemGll, countElemGll;
	startElemGll.push_back(temp_startElem);
	countElemGll.push_back(temp_countElem[XMPI::rank()]);
	startElemGll.push_back(0);
	countElemGll.push_back(nPntEdge);
	startElemGll.push_back(0);
	countElemGll.push_back(nPntEdge);
	
	std::vector<size_t> startElemNuGll, countElemNuGll;
	startElemNuGll.push_back(temp_startElemNu);
	countElemNuGll.push_back(temp_countElemNu[XMPI::rank()]);
	startElemNuGll.push_back(0);
	countElemNuGll.push_back(2);
	startElemNuGll.push_back(0);
	countElemNuGll.push_back(nPntEdge);
	startElemNuGll.push_back(0);
	countElemNuGll.push_back(nPntEdge);
	

	int temp_startElem_ani, temp_startElemNu_ani;
	std::vector<int> temp_countElem_ani(XMPI::nproc(),0);
	std::vector<int> temp_countElemNu_ani(XMPI::nproc(),0);

	XMPI::gather(elems_proc_ani, temp_countElem_ani, true); //non fourier fields
	temp_startElem_ani = std::accumulate(temp_countElem_ani.begin(), temp_countElem_ani.begin()+XMPI::rank(), 0);

	XMPI::gather(elemNus_proc_ani, temp_countElemNu_ani, true); //fourier fields 
	temp_startElemNu_ani = std::accumulate(temp_countElemNu_ani.begin(), temp_countElemNu_ani.begin()+XMPI::rank(),0);
	
	std::vector<size_t> startElem_ani, countElem_ani;
	startElem_ani.push_back(temp_startElem_ani);
	countElem_ani.push_back(temp_countElem_ani[XMPI::rank()]);
	
	std::vector<size_t> startElemGll_ani, countElemGll_ani;
	startElemGll_ani.push_back(temp_startElem_ani);
	countElemGll_ani.push_back(temp_countElem_ani[XMPI::rank()]);
	startElemGll_ani.push_back(0);
	countElemGll_ani.push_back(gll_ani);
	startElemGll_ani.push_back(0);
	countElemGll_ani.push_back(gll_ani);
	
	std::vector<size_t> startElemNuGll_ani, countElemNuGll_ani;
	startElemNuGll_ani.push_back(temp_startElemNu_ani);
	countElemNuGll_ani.push_back(temp_countElemNu_ani[XMPI::rank()]);
	startElemNuGll_ani.push_back(0);
	countElemNuGll_ani.push_back(2);
	startElemNuGll_ani.push_back(0);
	countElemNuGll_ani.push_back(gll_ani);
	startElemNuGll_ani.push_back(0);
	countElemNuGll_ani.push_back(gll_ani);

	// ---------- </DEFINE WRITING OFFSETS> ----------

	// ---------- <DEFINE NETCDF FILE> ----------
	NetCDF_Writer nc_writer = NetCDF_Writer();	
	if (XMPI::root()) { //Kuangdai implemented that creation and definition has to be serial
		
	std::vector<size_t> dimsElem, dimsElemGll, dimsElemNuGll;
	
	dimsElem.push_back(elems_all);	
	dimsElemGll.push_back(elems_all);
	dimsElemGll.push_back(nPntEdge);
	dimsElemGll.push_back(nPntEdge);	
	dimsElemNuGll.push_back(elemNus_all);
	dimsElemNuGll.push_back(2); // real/imag
	dimsElemNuGll.push_back(nPntEdge);
	dimsElemNuGll.push_back(nPntEdge);
	
	std::vector<size_t> dimsElem_ani, dimsElemGll_ani, dimsElemNuGll_ani;
	
	dimsElem_ani.push_back(elems_all_ani);	
	dimsElemGll_ani.push_back(elems_all_ani);
	dimsElemGll_ani.push_back(gll_ani);
	dimsElemGll_ani.push_back(gll_ani);
	dimsElemNuGll_ani.push_back(elemNus_all);
	dimsElemNuGll_ani.push_back(2);
	dimsElemNuGll_ani.push_back(gll_ani);
	dimsElemNuGll_ani.push_back(gll_ani);


	nc_writer.open(mFileName,true);	
	nc_writer.defModeOn();
	nc_writer.defineVariable<int>("Nus", dimsElem);
	nc_writer.defineVariable<int>("Nrs", dimsElem);
	nc_writer.defineVariable<int>("domain_decomposition", dimsElem);
	nc_writer.defineVariable<int>("sem_mesh", dimsElemGll);
	nc_writer.defineVariable<Real>("mesh_S", dimsElemGll);
	nc_writer.defineVariable<Real>("mesh_Z", dimsElemGll);
	nc_writer.defineVariable<Real>("mesh_vp", dimsElemNuGll);
	nc_writer.defineVariable<Real>("mesh_vs", dimsElemNuGll);

	nc_writer.defModeOff();
	
	nc_writer.writeVariableWhole("domain_decomposition", mMesh->mMsgInfo->mElemToProc);
	
	nc_writer.close();
	}
	// ---------- </DEFINE NETCDF FILE> ----------
	
	// ---------- <TEST ARRAY VS VECTOR> ----------	

	vec_ar2_RMatPP test_ar(elemNus_proc_ani,zero_ar2_RMatPP);
	vec_RMatPP zero_vec2(2, RMatPP::Zero());
	std::vector<vec_RMatPP> test_vec(elemNus_proc_ani, zero_vec2);

	for (int i=0;i<elemNus_proc_ani;i++)
		for (int j=0;j<2;j++)
			for (int ipol=0;ipol<=nPol;ipol++)
				for (int jpol=0;jpol<=nPol;jpol++) {
					test_ar[i][j](ipol,jpol)=1.;
					test_vec[i][j](ipol,jpol)=1.;
				}
	// ---------- </TEST ARRAY VS VECTOR> ----------	
				
	// ---------- <WRITE FIELDS> ----------	
	//  parallel NetCDF required
	nc_writer.openParallel(mFileName);

	nc_writer.writeVariableChunk("Nus", Nus, startElem, countElem);
	nc_writer.writeVariableChunk("Nrs", Nrs, startElem, countElem);
	nc_writer.writeVariableChunk("sem_mesh", mMesh->mMsgInfo->mLocElemToGlobPoints, startElemGll, countElemGll);
	nc_writer.writeVariableChunk("mesh_S", mesh_S, startElemGll, countElemGll);
	nc_writer.writeVariableChunk("mesh_Z", mesh_Z, startElemGll, countElemGll);
	nc_writer.writeVariableChunk("mesh_vp", mesh_vp, startElemNuGll, countElemNuGll);
	nc_writer.writeVariableChunk("mesh_vs", mesh_vs, startElemNuGll, countElemNuGll);


	nc_writer.close();
	// ---------- </WRITE FIELDS> ----------	


	
}