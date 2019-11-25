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
		for(int i = 0 ; i < 1000 ; i++){
			Add_Solid_Circle_H(i);
		}
	}

	void Sync_Simulation_And_Visualization_Data()
	{


		int count = 0;
		for (int i = 0; i < fluid.particles.Size(); i++)
		{
			OpenGLColor my_blue = OpenGLColor(0.0f, 0.0f, 1.f, 1.f);
			OpenGLColor my_red = OpenGLColor(1.0f, 0.0f, 0.f, 1.f);
			OpenGLColor my_yellow = OpenGLColor(1.0f, 1.0f, 0.f, 0.5f);
			auto opengl_circle = opengl_circles[i];
			opengl_circle->pos = V3(fluid.particles.X(i));
			if (fluid.particles.H(i) == 0)
			{
				opengl_circle->color = my_blue; // L region
			}
			else if (fluid.particles.H(i) == 1)
			{
				opengl_circle->color = my_yellow; // surface
				count++;
			}
			else if (fluid.particles.H(i) == 2)
			{
				opengl_circle->color = my_red; // boundary
				count++;
			}

			opengl_circle->Set_Data_Refreshed();
		}

		// std::cout <<"cnt:" << count <<std::endl;
		// std::cout << "parents:" << fluid.parents_idx.size() <<std::endl;
		// std::cout << fluid.new_particles_H.Size() << std::endl;
		for (int j = 0 ; j < 1000 ; j++){
			auto opengl_circle = opengl_circles_H[j];
			opengl_circle->visible = false;
		}
		for (int i = 0; i < fluid.new_particles_H.Size(); i++)
		{
			// OpenGLColor my_blue = OpenGLColor(0.0f, 0.0f, 1.f, 1.f);
			// OpenGLColor my_red = OpenGLColor(1.0f, 0.0f, 0.f, 1.f);
			OpenGLColor my_black = OpenGLColor(0.0f, 0.0f, 0.f, 1.f);
			auto opengl_circle = opengl_circles_H[i];
			opengl_circle->pos = V3(fluid.new_particles_H.X(i));
			opengl_circle->radius=fluid.new_particles_H.R(i) ;
			// std::cout << "r:" << fluid.new_particles_H.R(i)  << std::endl;
			// std::cout << "R2:" << fluid.particles.R(i)  << std::endl;
			opengl_circle->visible = true;
			
			opengl_circle->color = my_black; // surface
			
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
