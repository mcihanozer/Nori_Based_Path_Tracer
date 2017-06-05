/*
	Base class for inherit and create geometric objects such as
	Sphere.
*/

#if !defined(__GEOMETRICOBJECTS_H)
#define __GEOMETRICOBJECTS_H

#include <nori\mesh.h>
#include <nori\bsdf.h>

NORI_NAMESPACE_BEGIN


/*	
	DEVELOPER: MCO

	Base class for basic geometric objects (like sphere). it inhertis Mesh class
	of Nori and rearranges class methods. Override some of them.

	Using this class is essential for adding new classes wanted to render. Because,
	all rendered objects should be added into KD-tree of the scene. And if their
	base class is not Mesh, KD-tree could not put them in the tree and check for
	ray intersection
*/
class GeometricShape : public Mesh
{
public:
	
	//void samplePosition(const Point2f &sample, Point3f &p, Normal3f &n) const{	/*	NO DEFINITOIN FOR YOU!	*/	}

	//inline float pdf() const{	/*	NO DEFINITOIN FOR YOU!	*/	}

	EClassType getClassType() const { return EShape; }

	EMeshType getMeshType() const { return GeometricType; }

	inline const Transform &getTransform() const {	return m_transform;	}

protected:

	Transform	m_transform;	// Transformation matrix of the object
};

NORI_NAMESPACE_END

#endif