#ifndef __ImplicitGeometry_h__
#define __ImplicitGeometry_h__

template<int d> class ImplicitGeometry
{using VectorD=Vector<real,d>;
public:
	virtual real Phi(const VectorD& pos) const {return 0.;}
	virtual VectorD Normal(const VectorD& pos) const {return VectorD::Zero();}
	virtual VectorD get_center() const {return VectorD::Zero();}
	virtual real get_radius() const {return 1.0f;}
};

template<int d> class Bowl : public ImplicitGeometry<d>
{using VectorD=Vector<real,d>;
public:
	VectorD center;
	real radius;
	Bowl(VectorD _center=VectorD::Zero(),real _radius=1.):center(_center),radius(_radius){}
	virtual real Phi(const VectorD& pos) const {return radius-(pos-center).norm();}
	virtual VectorD Normal(const VectorD& pos) const {return (center-pos).normalized();}
	VectorD get_center() const {return center;}
	real get_radius() const {return radius;}
};

template<int d> class Sphere : public ImplicitGeometry<d>
{using VectorD=Vector<real,d>;
public:
	VectorD center;
	real radius;
	Sphere(VectorD _center=VectorD::Zero(),real _radius=1.):center(_center),radius(_radius){}
	virtual real Phi(const VectorD& pos) const {return (pos-center).norm()-radius;}
	virtual VectorD Normal(const VectorD& pos) const {return (pos-center).normalized();}
	VectorD get_center() const {return center;}
	real get_radius() const {return radius;}
};

#endif
