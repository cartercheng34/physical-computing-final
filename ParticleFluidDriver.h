//#####################################################################
// Particle Fluid Driver
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __ParticleFluidDriver_h__
#define __ParticleFluidDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "ParticleFluid.h"

template<int d> class ParticleFluidDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.02;
	ParticleFluid<d> fluid;
	Array<OpenGLSolidCircle*> opengl_circles;
	Array<OpenGLSolidCircle*> opengl_circles_H;

	Bowl<d>* bowl=nullptr;
public:
	virtual void Initialize()
	{
		////driver initialization, initialize simulation data
		real dx=.35;int nx=20;int ny=20;
		for(int i=0;i<nx;i++)
		{
			for(int j=0;j<ny;j++)
			{
				VectorD pos;
				pos[0]=(real)i*dx-1.;//initial positon
				pos[1]=(real)j*dx+3.;
				Add_Particle(pos);
		    }
		}

		// for(int i = 0 ; i < 700 ; i++){
		// 	VectorD pos;
		// 	pos[0]=0.0f;//initial positon
		// 	pos[1]=0.0f;
		// 	Add_Particle_H(pos, 1.0f/8, 1.0f/2);
		// }

		bowl=new Bowl<d>(VectorD::Unit(1)*8,8);
		fluid.env_objects.push_back(bowl);

		fluid.Initialize();
		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		if(bowl){
			auto opengl_circle=Add_Interactive_Object<OpenGLCircle>();
			opengl_circle->n=64;
			opengl_circle->pos=V3(bowl->center);
			opengl_circle->radius=bowl->radius;
			opengl_circle->color=OpenGLColor(1.f,.6f,.2f);
			opengl_circle->line_width=4.f;
			opengl_circle->Set_Data_Refreshed();
			opengl_circle->Initialize();}

		for(int i=0;i<fluid.particles.Size();i++){
			Add_Solid_Circle(i);}
		for(int i = 0 ; i < 700 ; i++){
			Add_Solid_Circle_H(i);
		}
	}

	void Sync_Simulation_And_Visualization_Data()
	{


		
		// auto record_idx = fluid.parents_idx;
		// auto delete_idx = fluid.prev_parents_idx;
		

		// for (int i = 0; i < fluid.parents_idx.size(); i++) {
		// 	bool find = false;
		// 	for (int j = 0; j < fluid.prev_parents_idx.size(); j++) {
		// 		if (fluid.parents_idx[i] == fluid.prev_parents_idx[j]) {
		// 			find = true;

		// 			break;
		// 		}

		// 	}
		// 	if (!find)
		// 		record_idx[i] = -1.0;
		// }

		// for (int i = 0; i < fluid.prev_parents_idx.size(); i++) {
		// 	bool find = false;
		// 	for (int j = 0; j < fluid.parents_idx.size(); j++) {
		// 		// std::cout << fluid.prev_parents_idx[j] << std::endl;
		// 		if (fluid.parents_idx[j] == fluid.prev_parents_idx[i]) {
		// 			find = true;

		// 			break;
		// 		}

		// 	}
		// 	if (!find)
		// 		delete_idx[i] = -1.0;
		// }

		// std::cout << "prev:" << fluid.particles_H.Size() << std::endl;
		// for (int i = 0 ; i < fluid.parents_idx.size() ; i++){
		// 	real L_radius = fluid.particles.R(fluid.parents_idx[i]);
		// 	real H_radius = L_radius / 8.0f;
		// 	real dx = 0.35f;
		// 	VectorD pos1 = fluid.particles.X(fluid.parents_idx[i]) + VectorD::Unit(0)*dx/4.0f + VectorD::Unit(1)*dx/4.0f;
		// 	VectorD pos2 = fluid.particles.X(fluid.parents_idx[i]) + VectorD::Unit(0)*dx/4.0f - VectorD::Unit(1)*dx/4.0f;
		// 	VectorD pos3 = fluid.particles.X(fluid.parents_idx[i]) - VectorD::Unit(0)*dx/4.0f + VectorD::Unit(1)*dx/4.0f;
		// 	VectorD pos4 = fluid.particles.X(fluid.parents_idx[i]) - VectorD::Unit(0)*dx/4.0f - VectorD::Unit(1)*dx/4.0f;
		// 	if (record_idx[i] == -1.0){
		// 		// add surface particles in fluid.particles_H
		// 		Add_Particle_H(fluid.particles.X(fluid.parents_idx[i]), fluid.particles.M(fluid.parents_idx[i]) , L_radius, fluid.parents_idx[i]);
		// 		// add small particles in fluid.new_particles_H
		// 		Add_Particle_small(pos1, fluid.particles.M(fluid.parents_idx[i]) / 8.0f, H_radius, fluid.parents_idx[i]);
		// 		Add_Particle_small(pos2, fluid.particles.M(fluid.parents_idx[i]) / 8.0f, H_radius, fluid.parents_idx[i]);
		// 		Add_Particle_small(pos3, fluid.particles.M(fluid.parents_idx[i]) / 8.0f, H_radius, fluid.parents_idx[i]);
		// 		Add_Particle_small(pos4, fluid.particles.M(fluid.parents_idx[i]) / 8.0f, H_radius, fluid.parents_idx[i]);

		// 		//draw physical particles
		// 		Add_Solid_Circle_H(i);
		// 		Add_Solid_Circle_H(i+1);
		// 		Add_Solid_Circle_H(i+2);
		// 		Add_Solid_Circle_H(i+3);
		// 	} // new particles that need to be drawn
				
			
		// }

		// opengl_circles_H.clear();
		// for (int i = 0 ; i < fluid.parents_idx.size() ; i++){

		// 	Add_Solid_Circle_H(i*4);
		// 	Add_Solid_Circle_H(i*4+1);
		// 	Add_Solid_Circle_H(i*4+2);
		// 	Add_Solid_Circle_H(i*4+3);
		// 	// opengl_circles_H[i*4]->Set_Data_Refreshed();
		// 	// opengl_circles_H[i*4+1]->Set_Data_Refreshed();
		// 	// opengl_circles_H[i*4+2]->Set_Data_Refreshed();
		// 	// opengl_circles_H[i*4+3]->Set_Data_Refreshed();
		// }
		// std::cout << fluid.parents_idx.size() << std::endl;
		// std::cout << opengl_circles_H.size() << std::endl;
		
		

		// std::cout << "prev:" << fluid.prev_parents_idx.size() << std::endl;
		// for (int i = fluid.prev_parents_idx.size() - 1; i >= 0; i--) {
			
		// 	if (delete_idx[i] == -1.0) {// old particles that need to be deleted
		// 		//delete surface
				
		// 		fluid.particles_H.X()->erase(fluid.particles_H.X()->begin() + i);
		// 		fluid.particles_H.V()->erase(fluid.particles_H.V()->begin() + i);
		// 		fluid.particles_H.F()->erase(fluid.particles_H.F()->begin() + i);
		// 		fluid.particles_H.C()->erase(fluid.particles_H.C()->begin() + i);
		// 		fluid.particles_H.R()->erase(fluid.particles_H.R()->begin() + i);
		// 		fluid.particles_H.P()->erase(fluid.particles_H.P()->begin() + i);
		// 		fluid.particles_H.D()->erase(fluid.particles_H.D()->begin() + i);
		// 		fluid.particles_H.I()->erase(fluid.particles_H.I()->begin() + i);
				
		// 		for(int j = 0 ; j < 4 ; j++){
		// 			// delete small 
		// 			fluid.new_particles_H.X()->erase(fluid.new_particles_H.X()->begin() + i*4);
		// 			fluid.new_particles_H.V()->erase(fluid.new_particles_H.V()->begin() + i*4);
		// 			fluid.new_particles_H.F()->erase(fluid.new_particles_H.F()->begin() + i*4);
		// 			fluid.new_particles_H.C()->erase(fluid.new_particles_H.C()->begin() + i*4);
		// 			fluid.new_particles_H.R()->erase(fluid.new_particles_H.R()->begin() + i*4);
		// 			fluid.new_particles_H.P()->erase(fluid.new_particles_H.P()->begin() + i*4);
		// 			fluid.new_particles_H.D()->erase(fluid.new_particles_H.D()->begin() + i*4);
		// 			fluid.new_particles_H.I()->erase(fluid.new_particles_H.I()->begin() + i*4);
		// 			opengl_circles_H.erase(opengl_circles_H.begin() + i*4);
		// 			// opengl_circles_H[i + j]->Set_Data_Refreshed();
		// 		}

		// 	}
		// }

		// std::cout << "deduct" << fluid.particles_H.Size() << std::endl;
		// std::cout << "surface:" << fluid.particles_H.Size() << std::endl;
		// std::cout << "small:" << fluid.new_particles_H.Size() << std::endl;


		for (int i = 0; i < fluid.particles.Size(); i++)
		{
			OpenGLColor my_blue = OpenGLColor(0.0f, 0.0f, 1.f, 1.f);
			OpenGLColor my_red = OpenGLColor(1.0f, 0.0f, 0.f, 1.f);
			OpenGLColor my_yellow = OpenGLColor(1.0f, 1.0f, 0.f, 1.f);
			auto opengl_circle = opengl_circles[i];
			opengl_circle->pos = V3(fluid.particles.X(i));
			if (fluid.particles.H(i) == 0)
			{
				opengl_circle->color = my_blue; // L region
			}
			else if (fluid.particles.H(i) == 1)
			{
				opengl_circle->color = my_yellow; // surface
			}
			else if (fluid.particles.H(i) == 2)
			{
				opengl_circle->color = my_red; // boundary
			}

			opengl_circle->Set_Data_Refreshed();
		}
		for (int j = 0 ; j < 700 ; j++){
			auto opengl_circle = opengl_circles_H[j];
			opengl_circle->visible = false;
		}
		for (int i = 0; i < fluid.new_particles_H.Size(); i++)
		{
			// OpenGLColor my_blue = OpenGLColor(0.0f, 0.0f, 1.f, 1.f);
			// OpenGLColor my_red = OpenGLColor(1.0f, 0.0f, 0.f, 1.f);
			OpenGLColor my_yellow = OpenGLColor(0.0f, 0.0f, 0.f, 1.f);
			auto opengl_circle = opengl_circles_H[i];
			opengl_circle->pos = V3(fluid.new_particles_H.X(i));
			opengl_circle->visible = true;
			
			opengl_circle->color = my_yellow; // surface
			
			opengl_circle->Set_Data_Refreshed();
		}
		
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		fluid.Advance(dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	////User interaction
	virtual bool Mouse_Click(int left,int right,int mid,int x,int y,int w,int h)
	{
		if(left!=1){return false;}
		Vector3f win_pos=opengl_window->Project(Vector3f::Zero());
		Vector3f pos=opengl_window->Unproject(Vector3f((float)x,(float)y,win_pos[2]));
		VectorD p_pos;for(int i=0;i<d;i++)p_pos[i]=(real)pos[i];
		real r=.1*static_cast<float>(rand()%1000)/1000.+.15;
		Add_Particle(p_pos);
		Add_Solid_Circle(fluid.particles.Size()-1);
		return true;
	}

protected:
	void Add_Particle(VectorD pos,real m=1.)
	{
		int i=fluid.particles.Add_Element();	////return the last element's index
		fluid.particles.X(i)=pos;
		fluid.particles.V(i)=VectorD::Zero();
		fluid.particles.R(i)=0.2f;
		fluid.particles.M(i)=m;
		fluid.particles.D(i)=1.;
		fluid.particles.C(i) = 0.5f;
		fluid.particles.I(i) = -1;
		fluid.particles.H(i) = 0;
	}

	void Add_Particle_H(VectorD pos, real m = 1., real radius = 1., int idx = 0, VectorD v = VectorD::Zero())
	{
		int i=fluid.new_particles_H.Add_Element();	////return the last element's index
		fluid.new_particles_H.X(i)=pos;
		fluid.new_particles_H.V(i)= v;
		fluid.new_particles_H.R(i)= radius;
		fluid.new_particles_H.M(i)=m;
		fluid.new_particles_H.D(i)=1.;
		fluid.new_particles_H.C(i) = 0.5f;
		fluid.new_particles_H.I(i) = idx; // record parent idx
		fluid.new_particles_H.H(i) = 0;
	}

	// void Add_Particle_small(VectorD pos, real m = 1., real radius = 1., int idx = 0, VectorD v = VectorD::Zero())
	// {
	// 	int i=fluid.new_particles_H.Add_Element();	////return the last element's index
	// 	fluid.new_particles_H.X(i)=pos;
	// 	fluid.new_particles_H.V(i)= v;
	// 	fluid.new_particles_H.R(i)= radius;
	// 	fluid.new_particles_H.M(i)=m;
	// 	fluid.new_particles_H.D(i)=1.;
	// 	fluid.new_particles_H.C(i) = 0.5f;
	// 	fluid.new_particles_H.I(i) = idx; // record parent idx
	// 	fluid.new_particles_H.H(i) = 0;
	// }

	void Add_Solid_Circle(const int i)
	{
		OpenGLColor c(0.5f,0.5f,0.5f,1.f);
		auto opengl_circle=Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circles.push_back(opengl_circle);
		opengl_circle->pos=V3(fluid.particles.X(i));
		opengl_circle->radius=fluid.particles.R(i);
		opengl_circle->color = OpenGLColor(fluid.particles.C(i), fluid.particles.C(i), fluid.particles.C(i),1.f);
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->Initialize();
	}

	void Add_Solid_Circle_H(const int i)
	{
		OpenGLColor c(0.0f,0.0f,0.0f,1.f);
		auto opengl_circle=Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circle->visible = false;
		opengl_circles_H.push_back(opengl_circle);
		opengl_circle->pos=V3(fluid.particles.X(i));
		opengl_circle->radius=fluid.particles.R(i) / 2.0;
		opengl_circle->color = c;
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->Initialize();
	}


	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif
